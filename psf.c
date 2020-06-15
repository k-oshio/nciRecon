//
//	psf of neural current field 
//	move this to RecKit framework later
//

#include "image.h"


void	test_1d();
void	test_2d_0();	//	sub-pixel psf
void	test_2d_1();	//	calc in k using analytical psf
void	test_2d_2();	//	point-by-point calc in k
void	test_3d();      //	3d psf
void	test_3d_2();    //	1d curve fitting
void	test_1d_rand();
void	test_hist();
void	test_white();	// white matter phantom

void	circle(KO_IMAGE	*f, float cy, float cx, float r, float a);
void	circle_var(KO_IMAGE	*f, float r, float a, float v);
void	rand_mg(float **p, float **q, int cy, int cx, int xflag, float a, float cor);
void	rand_ph(float **p, float **q, int cy, int cx, int xflag, float a, float cor);
void	bi_p(float **p, float **q, float cy, float cx, int xflag, float ph);
void	tri_m(float **p, float **q, float cy, float cx, int xflag, float ph);
void	tri_ph(float **p, float **q, float cy, float cx, int xflag, float ph);
void	add_1_r(KO_IMAGE *f, float cx, float cy, float a, int xflag);
void	add_lin(KO_IMAGE *f, int x, int y, float a, float ph, int xflag);
void	dipole(KO_IMAGE *f, float cy, float cx, float a, float th);

int
main(int ac, char *av[])
{
//	test_1d();
//	test_1d_rand();
//	test_2d_1();	// dipole phantom **** (current)
//	test_2d_2();	// 128 x 128, dft based
//	test_2d_0();	// PSF of single current dipole
//	test_3d();      // PSF of single current dipole, 3D psf integrated along slice
	test_3d_2();    // PSF of single current dipole, 1d curve fitting
//	test_hist();
//	test_white();

    return 0;
}

void
test_1d()
{
	int		i, j;
	float	a, x;
	float	re[1024];
	float	im[1024];
	int		n = 128;
	int		mode = 3;	// 0: sin(x), 1: tri-phase pulse, 2: step

	for (i = 0; i < n; i++) {
		re[i] = im[i] = 0;
		x = ((float)i - n/2) / n;
		switch (mode) {
		case 0 :
			a = sin(x * M_PI);
			if (i % 2 == 0) {
				a = -1.0;
			} else
			if (x > 0) {
				a = -a;
			}
			if (i == n/2) a = 0;
	//		a *= 1000;
			re[i] = a;
			im[i] = 0.0;
			break;
		case 1 :
			if (i == n/2) {
				re[i] = 1.0;
			}
			if (i == n/2 + 1) {
				re[i] = -0.5;
			}
			if (i == n/2 - 1) {
				re[i] = -0.5;
			}
			break;
		case 2 : // step
			if (i < n/2) {
				re[i] = 1.0;
			}
			if (i > n/2) {
				re[i] = -1.0;
			}
			if (i == n/2) {
				re[i] = 0.0;
			}
			break;
		case 3 :	// 1/r
			if (i < n/2) {
				re[i] = 1.0 / (n/2 - i);
			}
			if (i > n/2) {
				re[i] = -1.0 / (i - n/2);
			}
			if (i == n/2) {
				re[i] = 0.0;
			}
			break;
		case 4 :	// 1/r^2
			if (i < n/2) {
				re[i] = 1.0 / (n/2 - i)/(n/2 - i);
			}
			if (i > n/2) {
				re[i] = -1.0 / (i - n/2)/(i - n/2);
			}
			if (i == n/2) {
				re[i] = 0.0;
			}
			break;
		}
//printf("%d %f %f\n", i, re[i], im[i]);
	}
//exit(0);
	shift(re, im, n);
	fft(re, im, n, KO_FFT);
	shift(re, im, n);

	for (i = 0; i < n; i++) {
		printf("%d %f %f\n", i, re[i], im[i]);
	}
}

// analytical psf in k
void
test_2d_1()
{
	KO_IMAGE	*f;
	int			i, j, k, n = 128;	// index in k
	int			np = 3;
	float		x, y, a, th;
	float		**p, **q;
	float		nz = 0.5;

// constant part (image domain)
	f = new_image(n, n, KO_COMPLEX);
	circle(f, 1.0, 1.0, n * 0.37, 2.0);
	circle(f, 1.0, 1.0, n * 0.34, -1.2);

// FT
	shift2(f);
	fft2(f, KO_FFT);
	shift2(f);

// add dipoles (freq domain)
	a = 0.5;
	k = 0;
	n = np * 2 + 1;
	n *= n;
	for (i = -np; i <= np; i++) {
		for (j = -np; j <= np; j++, k++) {
			th = k * 2 * M_PI / n;
			dipole(f, i * 24/np, j * 24/np, a, th);
		}
	}

// add noise
	gwn(f, nz);

// k-space data
	put_image_block(&f, "PH_k.img", 1);

// FT
	shift2(f);
	fft2(f, KO_IFT);
	shift2(f);

// image domain data after 2D FT
	put_image_block(&f, "PH_im.img", 1);
}

// stochastic phantom
void
test_2d_2()
{
	KO_IMAGE	*f, *fn, *g;	// f:base, fn:cor_noise, g:k-data
	int			i, j, k, n = 128;	// index in k
	int			ii, jj;			// index in space
	float		x, y, d2, cs;
	float		**p, **q;
	float		**pp, **qq;
	float		wrx[256], wry[256];	// fourier kernel
	float		wix[256], wiy[256];	// fourier kernel
	float		re, im, wr, wi, a;
	int			mode = 0;	// 0:corr, 1:1/r

// constant part
	f = new_image(n, n, KO_COMPLEX);
	fn = new_image(n, n, KO_COMPLEX);
	circle(f, 0, 0, n * 0.37, 1.0);
	circle(f, 0, 0, n * 0.3, -0.6);

	g = new_image(n, n, KO_COMPLEX);
	pp = (float **)g->real;
	qq = (float **)g->imag;
	for (i = 0; i < n; i++) {
		printf("%d\r", i); fflush(stdout);
		for (k = 0; k < n; k++) {
			wry[k] = cos(M_PI * 2 * (i - n/2) * (k - n/2) / n);
			wiy[k] = sin(M_PI * 2 * (i - n/2) * (k - n/2) / n);
		}
		for (j = 0; j < n; j++) {
			// changes for each point in k
			cpy_image(f, fn);
			p = (float **)fn->real;
			q = (float **)fn->imag;
			switch (mode) {
			case 0 :	// correlation in space
				a = 1.0;
				rand_ph(p, q, n/2, n/2, 0, 0.1, 0.0);
				rand_ph(p, q, n/2 - 15, n/2 - 15, 1, a,  0.9);
				rand_ph(p, q, n/2 - 15, n/2 + 15, 0, a,  0.9);
				rand_ph(p, q, n/2 + 15, n/2 - 15, 1, a, -0.9);
				rand_ph(p, q, n/2 + 15, n/2 + 15, 0, a, -0.9);
				break;
			case 1 :	// 1/r (single pixel)
				break;
			}
			// calc 1pt fourier coef for (ki, kj)
			for (k = 0; k < n; k++) {
				wrx[k] = cos(M_PI * 2 * (j - n/2) * (k - n/2) / n);
				wix[k] = sin(M_PI * 2 * (j - n/2) * (k - n/2) / n);
			}
			re = im = 0;
			for (ii = 0; ii < n; ii++) {
				for (jj = 0; jj < n; jj++) {
					wr = wrx[jj] * wry[ii];
					wi = wix[jj] * wiy[ii];
					re += p[ii][jj] * wr - q[ii][jj] * wi;
					im += p[ii][jj] * wi + q[ii][jj] * wr;
				}
			}
			pp[i][j] = re;
			qq[i][j] = im;
		}
	}

	put_image_block(&f, "PH.img", 1);
	put_image_block(&g, "PH_kn.img", 1);
	shift2(g);
	fft2(g, KO_IFT);
	shift2(g);
	put_image_block(&g, "PH_ft.img", 1);
}

// single current dipole
void
test_2d_0()
{
	KO_IMAGE	*f;
	int			i, j, k, n = 512;
	float		x, y, r, cs, rrss, w;
	float		**p, **q;
	int			mode = 3;	// 0: forward, 1: inverse 2: PSF with filter 3: x / r
	float		sig = 0.2;
	float		th = M_PI * 0.1;

	f = new_image(n, n, KO_COMPLEX);

	p = (float **)f->real;
	q = (float **)f->imag;
	for (i = 0; i < n; i++) {
		y = i - n/2;
		for (j = 0; j < n; j++) {
			x = j - n/2;
			r = sqrt(x*x + y*y);
			if (r < 1.0) {
				p[i][j] = 0;
				q[i][j] = 0;
			} else {
				cs = x / r;
				switch (mode) {
				case 0 :
					p[i][j] = 0;
					q[i][j] = cs / (r*r);	// x / r^3
					break;
				case 1 :
					rrss = cs*cs + sig*sig;
					p[i][j] = 0;
					q[i][j] = cs / rrss;
					break;
				case 2 :
					rrss = cs*cs + sig*sig;
					p[i][j] = 0;
					q[i][j] = cos(atan2(y, x) + th) * cs / rrss;
					break;
				case 3 :
					p[i][j] = 0;
					q[i][j] = cs / r;	// x / r^2
					break;
				}
			}
		}
	}
	put_image_block(&f, "PSF_img.img", 1);
	shift2(f);
	fft2(f, KO_FFT);
	shift2(f);
	put_image_block(&f, "PSF_k.img", 1);
}

// single current dipole
void
test_3d()
{
	KO_IMAGE	*f;
	int			i, j, k, n = 512, nz = 16;   // 4x oversamping
	float		x, y, z, r, ph;
	float		**p, **q;

	f = new_image(n, n, KO_COMPLEX);

	p = (float **)f->real;
	q = (float **)f->imag;
	for (i = 0; i < n; i++) {
		y = i - n/2;
		for (j = 0; j < n; j++) {
			x = j - n/2;
            r = sqrt(x*x + y*y);
            if (r < 1.0) {
                p[i][j] = 0;
                q[i][j] = 0;
            } else {
                ph = 0;
                for (k = 0; k < nz; k++) {
                    z = k - nz/2;
                    r = sqrt(x*x + y*y + z*z);
                    ph += x / (r*r*r);
                }
                p[i][j] = 0;
                q[i][j] = ph;
            }
		}
	}
	put_image_block(&f, "PSF_img.img", 1);
	shift2(f);
	fft2(f, KO_FFT);
	shift2(f);
    crop(f);
    crop(f);    // factor of 4
	put_image_block(&f, "PSF_k.img", 1);

}

// 1D fit
void
test_3d_2()
{
    int     i, n = 128; //512;
    float   ph, x;

    for (i = 0; i < n; i++) {
        x = (float)i - n/2;
        if (x < 1) {
            ph = 0;
        } else {
            ph = 4000.0 * (0.5 * exp(-x / 10) + 0.5 * exp(-x / 50));
        }
        printf("%d %f\n", i, ph);
    }
}

void
test_hist()
{
	KO_IMAGE	*f;
	float		**p;
	int			i, j, n = 256;
	int			mode = 1;	// 0: orig, 1: single, 2: exp
	float		x, y, r0, r1;
	float		x0 = 1.0;	// sig center x
	float		y0 = 0;		// sig center y
	float		d0 = 0.5;	// 0
	float		d1 = 0.5;	// 1
	char		fname[256];

	f = new_image(n, n, KO_COMPLEX);
	p = (float **)f->real;
	for (i = 0; i < n; i++) {
		y = (n/2 - (float)i) * 4 / n;
		for (j = 0; j < n; j++) {
			x = ((float)j - n/2) * 4 / n;
			r0 = sqrt(x*x + y*y);
			r1 = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
			switch (mode) {
			case 0 :
				p[i][j] = (r0 - r1) / (r0 + r1);
				break;
			case 1 :
				p[i][j] = r0 / (r0 + r1);
				break;
			case 2 :
				p[i][j] = r0*r0 / (r0*r0 + r1*r1);
				break;
			case 3 :
				p[i][j] = (1 + exp(-r1 / d1)) * (1 - exp(-r0 / d0));
				break;
			}
		}
	}
	sprintf(fname, "HIST%d.img", mode);
	put_image_block(&f, fname, 1);
}

// white matter (space domain)
void
test_white()
{
	KO_IMAGE	*f;
	int			i, j, k, n = 1024;   // 4x oversamping
	float		a;
	float		**p, **q;

	f = new_image(n, n, KO_COMPLEX);

	p = (float **)f->real;
	q = (float **)f->imag;

	for (i = 0; i < 5; i++) {
		if (i % 2 == 0) {
			a = -2.0;
		} else {
			a = 3.0;
		}
		for (j = 0; j < 5; j++) {
			p[n/2 + i - 2][n/2 + j - 2] = a;
			q[n/2 + i - 2][n/2 + j - 2] = 0;
		}
	}

	put_image_block(&f, "PSF_img.img", 1);
	shift2(f);
	fft2(f, KO_FFT);
	shift2(f);
    crop(f);
    crop(f);    // factor of 4
	put_image_block(&f, "PSF_k.img", 1);
}


#define N 32

void
circle(KO_IMAGE	*f, float y, float x, float r, float a)
{
	float	**p;
	int		i, j, n;
	float	rr, xx, yy;

	p = (float **)f->real;
	n = f->xdim;
	r *= r;
	for (i = 0; i < n; i++) {
		yy = (i - n/2) - y;
		for (j = 0; j < n; j++) {
			xx = (j - n/2) - x;
			rr = xx * xx + yy * yy;
			if (rr < r) {
				p[i][j] += a;
			}
		}
	}
}

void
circle_var(KO_IMAGE	*f, float r, float a, float v)
{
	float	**p;
	int		i, j, n;
	float	rr, rx, ry;

	p = (float **)f->real;
	n = f->xdim;
	r *= r;
	for (i = 0; i < n; i++) {
		ry = (float)(i - n/2) * 20 * M_PI / n;
		for (j = 0; j < n; j++) {
			rx = (float)(j - n/2) * 20 * M_PI / n;
			rr = (i - n/2) * (i - n/2) + (j - n/2) * (j - n/2);
			if (rr < r) {
				p[i][j] += a + v * cos(rx) * cos(ry);
			}
		}
	}
}

void
bi_p(float **p, float **q, float cy, float cx, int xflag, float ph)
{
	int		i, ix, iy;
	float	r, a;
	float	cs, sn;

	ix = cx; iy = cy;
	for (i = -N; i < N; i++) {
		r = i - 0.5;
		a = ph / r;
		/* phase
		cs = cos(a);
		sn = sin(a);
		if (xflag) {
			p[iy    ][ix + i] = cs;
			p[iy + 1][ix + i] = cs;
			q[iy    ][ix + i] = sn;
			q[iy + 1][ix + i] = sn;
		} else {
			p[iy + i][ix    ] = cs;
			p[iy + i][ix + 1] = cs;
			q[iy + i][ix    ] = sn;
			q[iy + i][ix + 1] = sn;
		}
		*/
		// mag
		if (xflag) {
			p[iy    ][ix + i] += a;
			p[iy + 1][ix + i] += a;
		} else {
			p[iy + i][ix    ] += a;
			p[iy + i][ix + 1] += a;
		}
	}
}

// mag
void
tri_m(float **p, float **q, float cy, float cx, int xflag, float ph)
{
	int		i, ix, iy;
	float	r, a;
	float	cs, sn;

	ix = cx; iy = cy;
	for (i = -N; i < N; i++) {
		r = i;
		if (i == 0) {
			a = ph;
		} else {
			a = -0.5 * fabs(ph / r);
		}
//		cs = cos(a);
//		sn = sin(a);
		if (xflag) {
			p[iy    ][ix + i] += a;
			p[iy + 1][ix + i] += a;
		} else {
			p[iy + i][ix    ] += a;
			p[iy + i][ix + 1] += a;
		}
	}
}



// phase
void
tri_ph(float **p, float **q, float cy, float cx, int xflag, float ph)
{
	int		i, ix, iy;
	float	r, a;
	float	cs, sn;

	ix = cx; iy = cy;
	for (i = -N; i < N; i++) {
		r = i;
		if (i == 0) {
			a = ph;
		} else {
			a = -0.5 * fabs(ph / r);
		}
		cs = cos(a);
		sn = sin(a);
		if (xflag) {
			p[iy    ][ix + i] = cs;
			p[iy + 1][ix + i] = cs;
			q[iy    ][ix + i] = sn;
			q[iy + 1][ix + i] = sn;
		} else {
			p[iy + i][ix    ] = cs;
			p[iy + i][ix + 1] = cs;
			q[iy + i][ix    ] = sn;
			q[iy + i][ix + 1] = sn;
		}
	}
}

void
dipole(KO_IMAGE *f, float cy, float cx, float a, float th)
{
	int		i, j, n;
	float	**p, **q;
	float	cs, sn, tth, r, x, y, re, im;

	n = f->xdim;
	p = (float **)f->real;
	q = (float **)f->imag;
	
	for (i = 0; i < n; i++) {
		y = (i - n/2);
		for (j = 0; j < n; j++) {
			x = (j - n/2);
			r = x*x + y*y;
			r = sqrt(r);
			tth = atan2(y, x);
			cs = cos(tth + th);
	// psf of dipole at center
			re = cs * (0.5 * exp(-r / 10.0) + 0.5 * exp(-r / 50.0));
			im = 0;
	// add linear phase
			tth = (x*cx + y*cy) * 2 * M_PI / n;
			cs = cos(tth);
			sn = sin(tth);
			p[i][j] += a * (re * cs - im * sn);
			q[i][j] += a * (re * sn + im * cs);
		}
	}
}

KO_IMAGE *
auto_corr(KO_IMAGE *f)
{
	KO_IMAGE	*g;
	float		*p, *q;
	float		mg;
	int			i, j;
	int			xdim, ydim, size;

	xdim = f->xdim;
	ydim = f->ydim;
	size = f->size;
	g = new_image(xdim, ydim, KO_COMPLEX);
	cpy_image(f, g);
	shift2(g);
	fft2(g, KO_FFT);
	shift2(g);
	p = (float *)g->real[0];
	q = (float *)g->imag[0];
	for (i = 0; i < size; i++) {
		mg = p[i] * p[i] + q[i] + q[i];
		p[i] = mg;
		q[i] = 0;
	}
put_image_block(&g, "PSF_a_mg.img", 1);
	shift2(g);
	fft2(g, KO_IFT);	
	shift2(g);
	return g;
}

void
test_1d_rand()
{
	KO_IMAGE	*f, *g;
	float		**p, **q;
	float		a, c = 0.5;
	int			i, j;
	int			n = 128;

	f = new_image(n, n, KO_COMPLEX);

	p = (float **)f->real;
	q = (float **)f->imag;
	for (i = 0; i < n; i++) {
		a = nrml(1.0);
		p[i][0] = a;
		for (j = 1; j < n; j++) {
			a = nrml(1.0);
			p[i][j] = a * (1.0 - c) + p[i][j-1] * c;
		}
	}
	put_image_block(&f, "PSF_a.img", 1);
	
	g = auto_corr(f);
	put_image_block(&g, "PSF_a_ac.img", 1);

	shift2(f);
	fft2(f, KO_FFT);
	shift2(f);
	put_image_block(&f, "PSF_ft.img", 1);
	
}

void
rand_mg(float **p, float **q, int cy, int cx, int xflag, float a, float cor)
{
	int		i, j, w = 10;
	int		x, y;
	float	prev;

	cx -= w/2; cy -= w/2;

	if (xflag) {
		for (i = 0; i < w; i++) {
			j = 0;
			x = j + cx; y = i + cy;
			prev = nrml() * a;
			p[y][x] += prev;
			for (j = 1; j < w; j++) {
				x = j + cx; y = i + cy;
				p[y][x] += nrml() * a + prev * cor;
			}
		}
	} else {
		for (j = 0; j < w; j++) {
			i = 0;
			x = j + cx; y = i + cy;
			prev = nrml() * a;
			p[y][x] += prev;
			for (i = 1; i < w; i++) {
				x = j + cx; y = i + cy;
				p[y][x] += nrml() * a + prev * cor;
			}
		}
	}
}

void
rand_ph(float **p, float **q, int cy, int cx, int xflag, float a, float cor)
{
	int		i, j, w = 10;
	int		x, y;
	float	prev;

	cx -= w/2; cy -= w/2;

	if (xflag) {
		for (i = 0; i < w; i++) {
			j = 0;
			x = j + cx; y = i + cy;
			prev = nrml() * a;
			p[y][x] += prev;
			for (j = 1; j < w; j++) {
				x = j + cx; y = i + cy;
				q[y][x] += nrml() * a + prev * cor;
			}
		}
	} else {
		for (j = 0; j < w; j++) {
			i = 0;
			x = j + cx; y = i + cy;
			prev = nrml() * a;
			p[y][x] += prev;
			for (i = 1; i < w; i++) {
				x = j + cx; y = i + cy;
				q[y][x] += nrml() * a + prev * cor;
			}
		}
	}
}

// pos in space
// f is in freq domain
void
add_1_r(KO_IMAGE *f, float x, float y, float a, int xflag)
{
	int		i, j, jj, n;
	float	**p, **q;
	float	nz, ph, phy, phx;

	p = (float **)f->real;
	q = (float **)f->imag;
	n = f->xdim;

	if (xflag) {
		for (i = 0; i < n; i++) {
			phy = (float)y * i / n * 2 * M_PI;
			for (j = 0; j < n; j++) {
				nz = nrml();
				nz *= nz * a;
				ph = phy + (float)x * j / n * 2 * M_PI;
				if (j >= n/2) {
					p[i][j] +=  nz * cos(ph);
					q[i][j] +=  nz * sin(ph);
				}
				if (j < n/2) {
					p[i][j] += - nz * cos(ph);
					q[i][j] += - nz * sin(ph);
				}
			}
		}
	} else {
		for (j = 0; j < n; j++) {
			phx = (float)x * j / n * 2 * M_PI;
			for (i = 0; i < n; i++) {
				nz = nrml();
				nz *= nz * a;
				ph = phx + (float)y * i / n * 2 * M_PI;
				if (i >= n/2) {
					p[i][j] +=  nz * cos(ph);
					q[i][j] +=  nz * sin(ph);
				}
				if (i < n/2) {
					p[i][j] += - nz * cos(ph);
					q[i][j] += - nz * sin(ph);
				}
			}
		}
	}
}

float
sinc(float ph)
{
	float	a;
	if (ph == 0) {
		a = 1.0;
	} else {
		a = sin(ph) / ph;
	}
	return a;
}

float
one_r(float ph)
{
	float	a;
	if (ph == 0) {
		a = 0;
	} else {
		a = 1.0 / ph;
	}
	return a;
}

void
add_lin(KO_IMAGE *f, int x, int y, float a, float ph0, int xflag)
{
	int		i, j, n;
	float	**p, **q;
	float	ph;

	p = (float **)f->real;
	q = (float **)f->imag;
	n = f->xdim;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ph = 2.0 * M_PI * (j - n/2) / n;
		//	p[i][j] += a * (sinc(ph + ph0) - sinc(ph));	// X
			p[i][j] += a * sinc(ph + ph0);
		//	q[i][j] += a * one_r(ph + ph0);
//if (i == 0) {
//	printf("%d %f\n", j, p[i][j]);
//}
		}
	}
}




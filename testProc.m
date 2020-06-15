//
//	nciRec test routines
//

#import "RecKit.h"
#import "RecImageNCI.h"
#import "NumKit.h"

void		test1();	//
void		test2();	// t distribution
void		test3();	// F distribution
void		test4();	// PCA filter
void		test5();	// FT pair

RecImage    *mk_phantom(void);  // stripe
RecImage    *mk_phantom2(void); // normal [0, 1]
void		mk_phantom3(void);	// complex f-test

int
main()
{

//	mk_phantom3();
	
//	test1();
//	test2();
//	test3();
//	test4();	// PCA filter
	test5();	// FT pair

	return 0;
}

// === stat test ===
void
test1()
{
	RecImage	*stim, *ref;
//	RecImage	*m1, *m2;
	RecImage	*pimg;
	float		thres = 0.05;
	float		*p;
	int			i, n, cnt;
//	int			nbin = 100;

//	stim = [RecImage imageWithKOImage:@"phantom/PH_nz3.img"];
	stim = [RecImage imageWithKOImage:@"phantom/PH_nz2.img"];
	ref  = [RecImage imageWithKOImage:@"phantom/PH_nz1.img"];

//dump_hist(stim, 100, "hist_stim.dat");
//dump_hist(ref, 100, "hist_ref.dat");

	[ref copyLoopsOf:stim];

//[stim addConst:0.75];

	pimg = [stim pImageWithRef:ref complex:YES];


	[pimg saveAsKOImage:@"phantom/IMG_pimg0"];
	[pimg pvalThres:thres];
	[pimg saveAsKOImage:@"phantom/IMG_pimg1"];

	p = [pimg data];
	n = [pimg dataLength];
	cnt = 0;
	for (i = 0; i < n; i++) {
		if (p[i] > 0) {
			cnt++;
		}
	}
	printf("%f %f\n", thres, (float)cnt/n);

}

RecImage *
mk_phantom(void)
{
    RecImage    *img;
    int         i, j, k, l, n, m;
    int         x0, x, y;
    float       val;
    int         res = 1024;
    int         wd = 800;   // width of non-zero part
    int         hi = 600;   // height of non-zero part
    int         nblk = 10;
    int         blksz = wd / nblk;
    float       *p;

    img = [RecImage imageOfType:RECIMAGE_REAL xDim:res yDim:res];
    p = [img data];

    val = -1;   // background
    for (i = 0; i < hi; i++) {
        y = (res - hi) / 2 + i;
        for (j = 0; j < wd; j++) {
            x = (res - wd) / 2 + j;
            p[y * res + x] = val;
        }
    }
            
    for (k = 0; k < nblk; k++) {   // res block
        x0 = res/2 - wd/2 + k * blksz;   // left edge
        m = k + 1;  // width of single stripe
        for (i = 0; i < hi; i++) {
            n = blksz / m / 2;  // # of stripes
            for (j = 0; j < n; j++) {
                val = 1.0;
                for (l = 0; l < m; l++) {
                    x = x0 + j * m * 2 + l;
                    y = res/2 - hi/2 + i;
                    p[y * res + x] = val;
                }
            }
        }
    }


    return img;
}

RecImage *
mk_phantom2(void)
{
    RecImage    *img;
    int         i, j;
    float       val;
    int         res = 512;
    float       *p;
    RecNoise    *nz = [RecNoise noise];

    img = [RecImage imageOfType:RECIMAGE_REAL xDim:res yDim:res];
    p = [img data];

    for (i = 0; i < res; i++) {
        for (j = 0; j < res; j++) {
            val = [nz nrml];
            p[i * res + j] = val;
        }
    }

    return img;
}

void
mk_phantom3(void) // complex f-test
{
	RecImage	*img;
	int			i, j, k;
	float		val, r;
	int			res = 128;
	int			nAvg = 100;
	int			n;
	float		*p;
    RecNoise    *nz = [RecNoise noise];

    img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:res yDim:res zDim:nAvg];
    p = [img data];
	n = [img dataLength] * [img pixSize];
    for (i = 0; i < n; i++) {
		val = [nz nrml];	// mean = 0, sd = 1
		p[i] = val;
    }
	[img saveAsKOImage:@"phantom/PH_nz1.img"];
    for (i = 0; i < n; i++) {
		val = [nz nrml];
		p[i] = val;
    }
	[img saveAsKOImage:@"phantom/PH_nz2.img"];
    for (k = 0; k < nAvg; k++) {
    for (i = 0; i < res; i++) {
    for (j = 0; j < res; j++) {
		r = (i - res/2)*(i - res/2) + (j - res/2)*(j - res/2);
		r = sqrt(r);
		val = [nz nrml];
		if (r < 10.0) {
			val += 0.4;
		}
		p[k * res * res + i * res + j] = val;
    }
	}
	}
	[img saveAsKOImage:@"phantom/PH_nz3.img"];
}

void
test2()
{
	int		i, n = 100;
	float	x, t[3], nrml;
//float       Num_nrml(float m, float v);

// ==== t-dist
	for (i = 0; i < n; i++) {
		x = (float)(i - n/2) * 10 / n;
		t[0] = Num_tdist(x, 5);
		t[1] = Num_tdist(x, 10);
		t[2] = Num_tdist(x, 100);
		nrml = Num_gdist(x, 1.0);
		printf("%f %f %f %f %f\n", x, t[0], t[1], t[2], nrml);
	}
}

void
test3()
{
	int		i, n = 100;
	float	x, F[3];
	float	a = 400;
	float	b = 10;

// ==== F-dist
	for (i = 0; i < n; i++) {
		x = (float)i * 3 / n;
		F[0] = Num_fdist(x, 145, 145);
		printf("%f %f\n", x * b, F[0] * a);
	}
}

void
zeromean_vect(float *p, int n)
{
	float	mn;
	int		i;

	mn = 0;
	for (i = 0; i < n; i++) {
		mn += p[i];
	}
	mn /= n;
	for (i = 0; i < n; i++) {
		p[i] -= mn;
	}
}

void
normalize_vect(float *p, int n)
{
	float	mg;
	int		i;

	mg = 0;
	for (i = 0; i < n; i++) {
		mg += p[i] * p[i];
	}
	mg = 1.0 / sqrt(mg);

	for (i = 0; i < n; i++) {
		p[i] *= mg;
	}
}

void
dump_vect(float *p, int n)
{
	int		i;

	for (i = 0; i < n; i++) {
		printf("%d %f\n", i, p[i]);
	}
}

void
test4()
{
	RecImage	*img_U, *img_s, *img_E, *img_phs, *img_tmp;
	int			i, j, nAvg, nImg;
	float		*s, *U;
	float		mg, sg;
	float		*st1, *st2;
	float		nz, r1, r2;

	typedef struct {
		int		ix;		// original index
		float	s;		// sigma
		float	st1;	// stim1
		float	st2;	// stim2
	} sort_ent;
	int			(^compar)(const void *, const void *);
	sort_ent	*sortTab;

	img_U = [RecImage imageWithKOImage:@"images_sav/IMG_U"];
	img_s = [RecImage imageWithKOImage:@"images_sav/IMG_s"];
	img_E = [RecImage imageWithKOImage:@"images_sav/IMG_E"];
	img_phs = [RecImage imageWithKOImage:@"images_sav/IMG_time_0"];

// U freq analysis
	img_tmp = [img_U copy];
	[img_tmp fft1d:[img_U xLoop] direction:REC_INVERSE];
	[img_tmp saveAsKOImage:@"IMG_U_f.img"];

// k-filt E
	img_tmp = [img_E kFilt:0];
	[img_tmp saveAsKOImage:@"IMG_E_kFlt.img"];

	nAvg = [img_U xDim];
	nImg = [img_U yDim];
	U = [img_U data];
	s = [img_s data];

	st1 = (float *)malloc(sizeof(float) * nAvg);
	st2 = (float *)malloc(sizeof(float) * nAvg);

	sortTab = (sort_ent *)malloc(sizeof(sort_ent) * nAvg);
	
	// normalize s[]
//	normalize_vect(s, nAvg);

//dump_vect(s, nAvg);
	// normalize U[]
	for (i = 0; i < nImg; i++) {
		normalize_vect(U + i * nAvg, nAvg);
	}
//[img_U saveAsKOImage:@"IMG_U_n.img"];

	// Stim1
	for (i = 0; i < nAvg; i++) {
		switch (i % 3) {
		case 0 :	// ref
			st1[i] = 0; //-1;
			break;
		case 1 :	// st 1
			st1[i] = 1;
			break;
		case 2 :	// st 2
			st1[i] = 0;
			break;
		}
	}
	zeromean_vect(st1, nAvg);
	normalize_vect(st1, nAvg);
//dump_vect(st1, nAvg);

	// Stim2
	for (i = 0; i < nAvg; i++) {
		switch (i % 3) {
		case 0 :	// ref
			st2[i] = 0; //-1;
			break;
		case 1 :	// st 1
			st2[i] = 0;
			break;
		case 2 :	// st 2
			st2[i] = 1;
			break;
		}
	}
	zeromean_vect(st2, nAvg);
	normalize_vect(st2, nAvg);

	// U dot stim
	for (i = 0; i < nImg; i++) {
		mg = 0;
		for (j = 0; j < nAvg; j++) {
			sg = U[i * nAvg + j] * st1[j];
			mg += sg * sg;
		}
		sortTab[i].st1 = sqrt(mg);
	}
	for (i = 0; i < nImg; i++) {
		mg = 0;
		for (j = 0; j < nAvg; j++) {
			sg = U[i * nAvg + j] * st2[j];
			mg += sg * sg;
		}
		sortTab[i].st2 = sqrt(mg);
	}

	for (i = 0; i < nAvg; i++) {
		sortTab[i].ix = i;
		sortTab[i].s = s[i];
	}
	for (i = 0; i < nAvg; i++) {
//		printf("%d %f %f %f\n", sortTab[i].ix, sortTab[i].st1, sortTab[i].st2, s[i]);
	}

	// sort sftTab
	compar = ^int(const void *p1, const void *p2) {
		if ((((sort_ent *)p1)->st2 - ((sort_ent *)p1)->s) > (((sort_ent *)p2)->st2 - ((sort_ent *)p2)->s)) {
			return 1;
		} else 
		if ((((sort_ent *)p1)->st2 - ((sort_ent *)p1)->s) == (((sort_ent *)p2)->st2 - ((sort_ent *)p2)->s)) {
			return 0;
		} else {
			return -1;
		}
	};
	qsort_b(sortTab, nAvg, sizeof(sort_ent), compar);
	for (i = 0; i < nAvg; i++) {
//		printf("%d %d %f %f\n", i, sortTab[i].ix, sortTab[i].s - sortTab[i].st1, sortTab[i].s - sortTab[i].st2);
	}

// plot remaining noise energy
	nz = r1 = r2 = 1.0;
	for (i = 0; i < nAvg; i++) {
//		printf("%d %d %f %f %f %f %f\n", i, sortTab[i].ix, nz, r1, r2, r1 - nz, r2 - nz);
		nz -= (sortTab[i].s) * (sortTab[i].s);
		r1 -= (sortTab[i].st1) * (sortTab[i].st1);
		r2 -= (sortTab[i].st2) * (sortTab[i].st2);
	}

// PCA filter
	[img_E copyLoopsOf:img_phs];
	for (i = 0; i < 6; i++) {
		int ix = sortTab[i].ix;
		img_tmp = [img_E sliceAtIndex:ix];
		printf("s = %f\n", s[ix]);
		[img_tmp multByConst:s[ix] * 0.06];
	//	[img_tmp multByConst:32764.0];	// scale ????
		[img_phs addImage:img_tmp];
	}
[img_tmp saveAsKOImage:@"IMG_cpt"];
[img_phs saveAsKOImage:@"IMG_flt"];

	free(st1);
	free(st2);
	free(sortTab);
}

void
test5()
{
	RecImage	*tm, *fr;
	int			i;
	float		*p, *q, x;
	float		w = 1.0;
	int			fct = 4;
	int			nt = 64;
	int			nf = 1024;

	if (0) {	// ref
		tm = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nt];
		p = [tm real];
		for (i = 0; i < nt; i++) {
			x = (float)i;
			p[i] = (2 * x - x * x / w) * exp(-x / w);
		}

		fr = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nf];
		q = [fr real];
		for (i = 0; i < nt; i++) {
			q[i] = p[i];
		}
		[fr shift1d:[fr xLoop]];
		[fr saveAsKOImage:@"IMG_1"];
		
		[fr fft1d:[fr xLoop] direction:REC_INVERSE];
		[fr saveAsKOImage:@"IMG_2"];
	} else {	// test
		nt *= fct;
		w *= fct;
		tm = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nt];
		p = [tm real];
		for (i = 0; i < nt; i++) {
			x = (float)i;
			p[i] = (2 * x - x * x / w) * exp(-x / w);
		}

		nf *= fct;
		fr = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nf];
		q = [fr real];
		for (i = 0; i < nt; i++) {
			q[i] = p[i];
		}
		[fr shift1d:[fr xLoop]];
		[fr saveAsKOImage:@"IMG_1"];
		
		[fr fft1d:[fr xLoop] direction:REC_INVERSE];
		[fr saveAsKOImage:@"IMG_2"];
		[fr crop:[fr xLoop] to:[fr xDim]/fct];
		[fr saveAsKOImage:@"IMG_3"];
	}
}


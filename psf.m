//
//	current dipole imaging (CDI) numerical sim
//
//	=== plans ===
//	## clean-up ... 
//		pre-calculate cos/sin filter
//		phantom2() : direct dipole to k conversion (ft first, then k filter)
//	then...
//	iterative correction (with sparse image assumption)
//

#import "RecKit.h"
#import "NumKit.h"
#import "RecImageNCI.h"

// phantom
RecImage    *phantom(int dim);         // create phantom k data
void        circle(float *p, int xdim, float x, float y, float r, float a);
void        dipole(float *p, float *q, int dim, float cy, float cx, float a, float th);
void		quad(float *p, float *q, int dim);	// first ver. ###
void		make_ref(RecImage *f, float cy, float cx, float a, float th);
RecImage	*mkFilt(int dim, BOOL cosFlg);

// original k-filter
RecImage    *proc1(RecImage *img);  // extract imag part (phase)
RecImage    *kFilt(RecImage *img, int mode);  // k-space filter (cos/sin), k to dipole
void		make_sparse(RecImage *x, RecImage *y, float thres);
RecImage	*re_project(RecImage *x, RecImage *y);

// p-space / k-energy based methods
RecImage    *kShift(RecImage *img, int xsft, int ysft, BOOL flt); // shift in k
RecImage    *pSpectPt(RecImage *img, int cx, int cy, BOOL flt);	// calc p spectrum at (cx, cy)
RecImage    *pSpect(RecImage *img);								// calc p spectrum
RecImage    *pCirc(RecImage *img, int nth, BOOL flt);	// circular sampling in p-space
RecImage    *kCirc(RecImage *img, int nth, BOOL flt);	// circular sampling in k-space
RecImage    *pFilt(RecImage *img, int freq);

// k-space energy (side)
void        kSpEn(RecImage *img);

// tools
float       snr(RecImage *img);
void        dumpSig(RecImage *img1, RecImage *img2);
void		calc_cnr_ft(RecImage *f);
void		calc_cnr_kflt(RecImage *g);
void		calc_angle(RecImage *g);

void		test_quad();

int         mode = 0;       // 0: k-filter, 1:shift, 2:spect, 3:p-circ, p-filter
							// 8: iterative
float		mtf = 1.0; //0.5;		// phantom gen
float		rec_mtf = 2.0;		// recon filter
float       nz = 0.000; //0.001; //0.0001; //0.001; //0.004;//0.0008;    // 0.0008; //0.001;
float		rng = 0.9;

int
main()
{
    @autoreleasepool {
		int			i, n;
		float		rms, orms;
        RecImage    *ph;		// phantom data(k), read only
		RecImage	*ref;
		RecImage	*f, *g, *h;	// dipole(k)
		RecImage	*fx0, *fy0;	// dipole(k) initial estimate
		RecImage	*fx, *fy;	// dipole(k)
		RecImage	*cs, *sn;	// cos filter
		RecImage	*d;			// null space image
        BOOL        flt = YES;
		NSString	*path;

		printf("mode = %d\n", mode);

		test_quad();
		exit(0);

        system("rm PH_*");
        ph = phantom(128);             // k-space
        [ph saveAsKOImage:@"PH_k"];
		cs = mkFilt(128, YES);
		sn = mkFilt(128, NO);
//		ref = [RecImage imageWithKOImage:@"PH_ref"];	// created by phantom()
//		[ref removePointLoops];

        g = [ph copy];
        [g fft2d:REC_FORWARD];
		[g gauss2DLP:rec_mtf];	
        [g saveAsKOImage:@"PH_ft"];
	calc_cnr_ft(g);

		g = [ph copy];
		[g fft2d:REC_FORWARD];	// -> img
		[g removeRealPart];
//		[g saveAsKOImage:@"PH_imag"];

        switch (mode) {
        case 0 :    // k-filter (same as nciRec2)
			g = [g kFilt:0];
			[g gauss2DLP:rec_mtf];	
            [g saveAsKOImage:@"PH_kflt"];
		calc_cnr_kflt(g);
		calc_angle(g);
			[ph fft2d:REC_FORWARD];
			f = [ph fusionWith:g gain:2.0 mode:0];
			[f saveAsKOImage:@"PH_fus1"];
			f = [ph fusionWith:g gain:2.0 mode:1];
			[f saveAsKOImage:@"PH_fus2"];
            break;
        case 1 :    // null space
			g = [ph copy];
			[g fft2d:REC_FORWARD];	// -> img
			[g removeRealPart];
            [g saveAsKOImage:@"PH_imag"];

		//	g = [g kFilt];

            [g saveAsKOImage:@"PH_kflt"];
			[ph fft2d:REC_FORWARD];
			f = [ph fusionWith:g gain:2.0 mode:0];
			[f saveAsKOImage:@"PH_fus"];
            break;
        case 2 :    // p-spect
            g = pSpectPt(f, 88, 40, YES);   // peak
            [g saveAsKOImage:@"PH_spct_s"];
            g = pSpectPt(f, 88, 41, YES);   // side-lobe
            [g saveAsKOImage:@"PH_spct_ss"];
            g = pSpectPt(f, 83, 46, YES);     // noise
            [g saveAsKOImage:@"PH_spct_n"];
			g = pSpect(f);
			[g saveAsKOImage:@"PH_spct"];
            break;
        case 3 :    // correlation
        // 1st method
        //  g = kFilt(f, 0);                           // k-filt (2pt)
        //  g = pCirc(f, 32, YES);    // p-filt
            g = kCirc(f, 48, YES); [g saveAsKOImage:@"PH_pc"];
            g = pFilt(g, 1); [g saveAsKOImage:@"PH_pf"];    // ke-filt (32 pt)
 
        // 2nd method
            h = kFilt(f, 0);                               // k-filt (2pt)
        //  h = pCirc(f, 32, YES); h = pFilt(h, 1);    // p-filt
        //    h = kCirc(f, 2, YES); h = pFilt(h, 1);    // ke-filt (32 pt)
            [h saveAsKOImage:@"PH_kf"];

            [g scaleToVal:1.0];
            [h scaleToVal:1.0];

            printf("p-Filt:%6.4f\n", snr(g));           // p-Filt:5.08, pk-Filt:5.83
            printf("k-Filt:%6.4f\n", snr(h));           // k-Filt:5.40
            dumpSig(h, g);

            [h addImage:g];
            printf("pk:%6.4f\n", snr(h));               // pk:7.23(7.41),     ppk:7.75(7.95)
            [h saveAsKOImage:@"PH_pk"];

            break;
        case 4 :    // circular sampling in p-space
            g = kCirc(f, 32, YES);
            [g saveAsKOImage:@"PH_kcirc"];
            g = pFilt(g, 1);
            [g saveAsKOImage:@"PH_kpflt"];
            printf("kp-Filt:%6.4f\n", snr(g));          // kp-Filt:5.12
        case 5 :    // k-space energy (p-grad)
            f = [RecImage imageWithKOImage:@"../gr.raw"];
            [f removePointLoops];
            kSpEn(f);
            break;
        case 6 :    // real data test
            f = [RecImage imageWithKOImage:@"../testavg.raw"];
			f = proc1(f);
		//	f = [f sliceAtIndex:2];
        //    g = pSpectPt(f, 30, 64, YES);
            g = pSpect(f);
            [g saveAsKOImage:@"IMG_spct"];
            break;
        case 7 :    // real data test
            f = [RecImage imageWithKOImage:@"../testavg.raw"];
			f = proc1(f);
			// kfilt / pfilt
            g = kCirc(f, 2, YES);	// 16
            [g saveAsKOImage:@"IMG_kcirc"];
            g = pFilt(g, 1);
            [g saveAsKOImage:@"IMG_pfilt"];
            break;
        case 8 :    // iterative
		// input is image data (imaginary only)
			f = [g copy];
			[f fft2d:REC_INVERSE];
[f fGauss2DLP:rec_mtf];	
			// cos/sin filter (initial est)
            fx = [f copy]; [fx scaleByImage:cs]; [fx saveAsKOImage:@"PH_kx"];
			fy = [f copy]; [fy scaleByImage:sn]; [fy saveAsKOImage:@"PH_ky"];

			// save f0
			fx0 = [fx copy];
			fy0 = [fy copy];

			//	FT
			[fx fft2d:REC_FORWARD]; [fx saveAsKOImage:@"PH_x"];	// img
			[fy fft2d:REC_FORWARD]; [fy saveAsKOImage:@"PH_y"];
			g = [fx copy];
			[g copyRealOf:fy];
			[g saveAsKOImage:@"PH_xy"];

			n = 10;
			rms = 0;
			for (i = 0; i < n; i++) {
			@autoreleasepool {
			//	printf("iter = %d\n", i);
				// make sparse
			//	make_sparse(fx, fy, (float)i * 0.7 + 1.0);
				make_sparse(fx, fy, 3.0);

				[fx removeImagPart];
				[fy removeImagPart];
[fx fGauss2DLP:rec_mtf];	
[fy fGauss2DLP:rec_mtf];	
			[fx saveAsKOImage:@"PH_x0"];
			[fy saveAsKOImage:@"PH_y0"];

				[fx fft2d:REC_INVERSE];		// F''x
				[fy fft2d:REC_INVERSE];		// F''y
				[fx subImage:fx0];			// (F''x - F'x)
				[fy subImage:fy0];			// (F''y - F'y)
			[fx saveAsKOImage:@"PH_kx0"];
			[fy saveAsKOImage:@"PH_ky0"];
				[fx scaleByImage:sn];
				[fy scaleByImage:cs];
				h = [fx copy];
				[h subImage:fy];			// d
				[h saveAsKOImage:@"PH_reproj"];
			orms = rms;
			rms = [h rmsVal];
			printf("%d %e\n", i, rms - orms);
				// correct
				fx = [h copy]; [fx scaleByImage:sn];
				fy = [h copy]; [fy scaleByImage:cs];
				[fx addImage:fx0];
				[fy subImage:fy0];
				[fy negate];
			[fx saveAsKOImage:@"PH_kx1"];
			[fy saveAsKOImage:@"PH_ky1"];
				[fx fft2d:REC_FORWARD]; 
				[fy fft2d:REC_FORWARD]; 
				g = [fx copy];
				[g copyRealOf:fy];
				path = [NSString stringWithFormat:@"PH_xy%03d", i + 1];
				[g saveAsKOImage:path];
			// rms error
			//	[g copyLoopsOf:ref];
			//	[g subImage:ref];
			//	rms = [g rmsVal];
			//	printf("%d %f\n", i, rms);
		//	[g saveAsKOImage:@"PH_error"];
		//	exit(0);
			}
			}
            break;
        }
        return 0;
    }   // autoreleasepool
}

RecImage *
proc1(RecImage *f)
{
    RecImage    *img = [f copy];
    [img fft2d:REC_FORWARD];
    [img removeRealPart];
    [img fft2d:REC_INVERSE];
    return img;
}

RecImage *
kFilt(RecImage *f, int mode)
{
    RecImage    *x, *y;
    x = [f copy];
    [x sinFilt:mode ang: 0];
    [x fft2d:REC_FORWARD];
	y = [f copy];
    [y sinFilt:mode ang: M_PI/2.0];
    [y fft2d:REC_FORWARD];
    [x copyRealOf:y];
    return x;
}

void
make_sparse(RecImage *x, RecImage *y, float thres)
{
	float	*xp = [x data];
	float	*yp = [y data];
	float	mg, mn;
	int		i, len = [x xDim] * [x yDim];

// find mean
	mn = 0;
	for (i = 0; i < len; i++) {
		mg = xp[i] * xp[i] + yp[i] * yp[i];
		mn += mg;
	}
	mn /= len;
	mn = sqrt(mn);
	thres *= mn;
	for (i = 0; i < len; i++) {
		mg = xp[i] * xp[i] + yp[i] * yp[i];
		if (sqrt(mg) < thres) {
			xp[i] = 0;
			yp[i] = 0;
		}
	}
}

RecImage *
re_project(RecImage *x, RecImage *y)
{
	RecImage	*img;

	return img;
}

// p-space encoding
RecImage *
kShift(RecImage *img, int xsft, int ysft, BOOL flt)
{
    RecImage		*sft = [img copy];
	RecLoopControl	*lc = [sft control];
    int				xDim = [img xDim];
    int				yDim = [img yDim];
    int				dataLength = [img dataLength];
    float			*p, *q;
	int				nSlc;
    int				i, j, k, x, y, ix;
    float			w, sd;

    sd = xDim * 0.1;

	[lc rewind];
	[lc deactivateXY];
	nSlc = [lc loopLength];

	for (k = 0; k < nSlc; k++) {
		p = [sft currentDataWithControl:lc];
		q = p + [sft dataLength];
		for (i = ix = 0; i < yDim; i++) {
			for (j = 0; j < xDim; j++, ix++) {
				if (flt) {
				// gaussian window
					x = j - xDim/2 - xsft;
					y = i - yDim/2 - ysft;
					w = x*x + y*y;
					w = exp(-w/(2.0*sd*sd));
					p[ix] *= w;
					q[ix] *= w;
				} else {
					x = j - xsft;
					y = i - ysft;
					if (y < 0 || y >= yDim || x < 0 || x >= xDim) {
						p[ix] = 0;
						q[ix] = 0;
					}
				}
			}
		}
		[lc increment];
	}
    return sft;
}

RecImage *
pSpectPt(RecImage *img, int cx, int cy, BOOL flt)     // calc p spectrum at (cx, cy)
{
    RecImage    *psp, *tmp;
    float       *p, *q;
    float       *pp, *qq;
    int         i, j, xDim, yDim;
    int         xsft, ysft, ofs;
    int         dim = 32;

    xDim = [img xDim];  // dim of src
    yDim = [img yDim];  // dim of src

    ofs = cy * xDim + cx;

    psp = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:dim yDim:dim];
    p = [psp data];
    q = p + [psp dataLength];

    for (i = 0; i < dim; i++) {
        @autoreleasepool {
            ysft = (float)i*yDim/dim - yDim/2;
            for (j = 0; j < dim; j++) {
                xsft = (float)j*xDim/dim - xDim/2;
                tmp = kShift(img, xsft, ysft, flt);
                [tmp fft2d:REC_FORWARD];
                pp = [tmp data] + ofs;
                qq = pp + [tmp dataLength];
                p[i*dim + j] = *pp;
                q[i*dim + j] = *qq;
            }
        }
    }

    return psp;
}

RecImage *
pSpect(RecImage *img)     // calc p spectrum
{
    RecImage		*psp, *tmp;
	RecLoopControl	*lc;
	NSMutableArray	*loops;
	RecLoop			*pxLp, *pyLp;
	RecLoopIndex	*pxLi, *pyLi;
    float			*p, *q;
    float			*pp, *qq;
    int				i, j, xDim, yDim;
    int				xsft, ysft, ofs;
    int				dim = 32;
	BOOL			flt = YES;

    xDim = [img xDim];  // dim of src
    yDim = [img yDim];  // dim of src

	loops = [NSMutableArray arrayWithArray:[img loops]];
	pxLp = [RecLoop loopWithDataLength:dim];
	pyLp = [RecLoop loopWithDataLength:dim];
	[loops addObject:pyLp];
	[loops addObject:pxLp];
	psp = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loops];
	lc = [psp control];
	[lc deactivateXY];
	pxLi = [lc loopIndexForLoop:pxLp];
	pyLi = [lc loopIndexForLoop:pyLp];

	// py
    for (i = 0; i < dim; i++) {
        @autoreleasepool {
            ysft = (float)i*yDim/dim - yDim/2;
			[pyLi setCurrent:i];
			// px
            for (j = 0; j < dim; j++) {
                xsft = (float)j*xDim/dim - xDim/2;
				[pxLi setCurrent:j];
				// tmp: [slc, y, x]
                tmp = kShift(img, xsft, ysft, flt);
                [tmp fft2d:REC_FORWARD];
				// copy tmp to psp
				[psp copyImage:tmp withControl:lc];
            }
        }
    }

    return psp;
}

RecImage *
pCirc(RecImage *img, int nth, BOOL flt)             // circular sampling in p-space
{
    RecImage        *crc, *tmp;
    RecLoop         *ang;
    RecLoopIndex    *li;
    RecLoopControl  *lc;
    NSMutableArray  *loopArray;
    int             i, xDim = [img xDim];
    float           th, r;
    int             x, y;

    ang = [RecLoop loopWithDataLength:nth];
    loopArray = [NSMutableArray arrayWithArray:[img loops]];
    [loopArray insertObject:ang atIndex:[loopArray count] - 2];
    crc = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    lc = [crc control];
    li = [lc loopIndexForLoop:ang];
    [lc deactivateLoop:ang];

    r = [img xDim] * 0.5;
    for (i = 0; i < nth; i++) {
        th = (float)i * 2 * M_PI / nth;
        x = cos(th) * r;
        y = sin(th) * r;
    //    printf("%f %d %d\n", th, x, y);
        tmp = kShift(img, x, y, flt);
        [tmp fft2d:REC_FORWARD];
        [li setCurrent:i];
        [crc copyImage:tmp withControl:lc];
    }

    return crc;
}

RecImage *
kCirc(RecImage *img, int nth, BOOL sym)             // circular sampling in p-space
{
    RecImage        *crc, *tmp;
    RecLoop         *ang;
    RecLoopIndex    *liAng;
    RecLoopControl  *lc, *lcSlc;
	int				loopLen;
    NSMutableArray  *loopArray;
    int             i, j, xDim = [img xDim], yDim = [img yDim];
    int             dataLength = [img dataLength];
    float           th, th2, *p, *q, w;
    int             x, y, sign;
    int             mode = 0;   // 0: cos, 1: rect

    if (nth < 0) {
        sign = -1;
        nth = -nth;
    } else {
        sign = 1;
    }
    ang = [RecLoop loopWithDataLength:nth];
    loopArray = [NSMutableArray arrayWithArray:[img loops]];
    [loopArray insertObject:ang atIndex:[loopArray count] - 2];
    crc = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    lc = [crc control];
	lcSlc = [RecLoopControl controlWithControl:lc];
    liAng = [lc loopIndexForLoop:ang];
    [lc deactivateLoop:ang];
	[lc deactivateXY];
	[lcSlc deactivateAll];
	[lcSlc activateXY];
	loopLen = [lc loopLength];
	tmp = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[img yLoop], [img xLoop], nil];
	p = [tmp real];
	q = [tmp imag];

    for (i = 0; i < nth; i++) {
        th = sign * (float)i * M_PI / nth;  // was 2_PI
		[liAng setCurrent:i];
		for (j = 0; j < loopLen; j++) {	// outerloop (slice)
			[tmp copyImage:img withControl:lcSlc];
			for (y = 0; y < yDim; y++) {
				for (x = 0; x < xDim; x++) {
					th2 = atan2((float)y - yDim/2, (float)x - xDim/2);
					w = cos(th2 - th);
					if (mode == 1) {
						if (w > 0) {
							w = 1.0;
						} else {
							w = -1.0;
						}
					}
					p[y*xDim + x] *= w; // better than rect (w = 1.0)
					q[y*xDim + x] *= w;
					if (!sym) {
						if (w < 0) {
							p[y*xDim + x] = 0;
							q[y*xDim + x] = 0;
						}
					}
				}
			}
			[tmp fft2d:REC_FORWARD];
			[crc copyImage:tmp withControl:lcSlc];

			[lc increment];
		}
    }

    return crc;
}

// multi-slice doesn't work ####
RecImage *
pFilt(RecImage *crc, int freq)
{
    RecImage        *flt;   // result
    RecLoop         *ang = [crc loopAtIndex:[crc dim] - 3];
	NSMutableArray	*loopArray;
    RecLoopControl  *lc = [crc control];
    int             dataLength;
    int             i, j, nth = [ang dataLength];
    int             skip = [crc skipSizeForLoop:ang];
    float           *buf = (float *)malloc(sizeof(float) * nth);
    float           *p, *q, *pp;
    float           cs, sn, th;

	loopArray = [NSMutableArray arrayWithArray:[crc loops]];
	[loopArray removeObject:ang];
	flt = [RecImage imageOfType:RECIMAGE_COMPLEX withLoopArray:loopArray];
    dataLength = [flt dataLength];
    p = [flt data];
    q = p + dataLength;

	[lc deactivateLoop:ang];
    [lc rewind];
    for (i = 0; i < dataLength; i++) {
        pp = [crc currentDataWithControl:lc];
        // copy ang data to buf
        for (j = 0; j < nth; j++) {
            buf[j] = *pp;
            pp += skip;
        }
        // project to sin/cos
        cs = sn = 0;
        for (j = 0; j < nth; j++) {
            th = freq * (float)j * M_PI / nth;  // was 2_PI
            cs += cos(th) * buf[j];
            sn += sin(th) * buf[j];
        }
        p[i] = cs;
        q[i] = sn;

        [lc increment];
    }
    free(buf);

    return flt;
}

RecImage *
phantom(int dim)
{
	RecImage	*f, *fxy;
    RecImage    *frot;
	float		*p, *q, th;
	int			i, j, k, len, n, np = 2;
    float       x, y;
	float		a = 1.0;
    int         mode = 2;   // 0: flat, 1: point, 2: grid

	f = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:dim yDim:dim];
	len = [f dataLength];
	p = [f data];
	q = p + len;

    [f setConst:1.0];
    [f sinFilt:0 ang:0.0];
//    [f saveAsKOImage:@"PH_filt"];

    [f clear];

// ==== space elements =====
	circle(p, dim, 1.0, 1.0, dim * 0.37,  2.0);
	circle(p, dim, 1.0, 1.0, dim * 0.34, -1.2);

// ==== freq elements ====
    [f fft2d:REC_INVERSE];

    switch (mode) {
    case 0 :    // flat
        break;
    case 1 :    // point
        x = y = 10;
        th = 0.0 * M_PI/2;
        dipole(p, q, dim, y, x, a, th);
        break;
    case 2 :    // grid
        k = 0;
        n = np * 2 + 1;
		n = n * n;
        for (i = -np; i <= np; i++) {
            for (j = -np; j <= np; j++, k++) {
                th = ((float)k - n/2) * 2 * M_PI * rng / n;
                x = (float)j * 24/np + 0.0;
                y = (float)i * 24/np + 0.0;
                dipole(p, q, dim, y, x, a, th);
            }
        }
    // UL marker
    //    x = (-np + 0.5) * 24/np;
    //    y = (-np + 0.5) * 24/np;
    //    th = 0;
    //    dipole(p, q, dim, y, x, a, th);
        
        break;
    }

	// gaussian filter (psf)
	[f fGauss2DLP:mtf];

    // add noise
    [f addGWN:nz relative:YES];

    return f;
}

void
circle(float *p, int dim, float x, float y, float r, float a)
{
	int		i, j, ix;
	float	rr, xx, yy;

	r *= r;
	for (i = 0; i < dim; i++) {
		yy = (i - dim/2) - y;
		for (j = 0; j < dim; j++) {
			xx = (j - dim/2) - x;
			rr = xx * xx + yy * yy;
			if (rr < r) {
                ix = i * dim + j;
				p[ix] += a;
			}
		}
	}
}

void
dipole(float *p, float *q, int n, float cy, float cx, float a, float th)
{
	int		i, j, ix;
	float	cs, sn, tth, r, x, y, re, im;

    a /= (n * n);   // after forward FT

	for (i = 0; i < n; i++) {
		y = (i - n/2);
		for (j = 0; j < n; j++) {
			x = (j - n/2);
			r = x*x + y*y;
			r = sqrt(r);
			tth = atan2(y, x);
			sn = sin(tth + th);
	// psf of diple at center
		//	re = cs * (0.5 * exp(-r / 10.0) + 0.5 * exp(-r / 50.0));
			re = sn;
			im = 0;
	// add linear phase
			tth = (x*cx + y*cy) * 2 * M_PI / n;
			cs = cos(tth);
			sn = sin(tth);
            ix = i * n + j;
			p[ix] += a * (re * cs - im * sn);
			q[ix] += a * (re * sn + im * cs);
		}
	}
}

void
quad(float *p, float *q, int dim)	// first ver. ###
{
	float		x, y;
	float		r, th;
	float		cs, sn;
	int			i, j;

	for (i = 0; i < dim; i++) {
		y = ((float)i - dim/2) / dim;
		for (j = 0; j < dim; j++) {
			x = ((float)j - dim/2) / dim;
			r = sqrt(x*x + y*y);
			th = atan2(y, x);
			if (r == 0) {
				p[i * dim + j] = 0;
			} else {
				p[i * dim + j] = cos(th * 2) / (r*r*r);
			}
		}
	}
}

void
test_quad()
{
	RecImage	*img;
	float		*p, *q;
	img = [RecImage imageOfType:RECIMAGE_REAL xDim:128 yDim:128];
	p = [img data];
	q = p + [img dataLength];

	quad(p, q, [img xDim]);
	[img saveAsKOImage:@"test_psf_quad.img"];
	[img fft2d:REC_INVERSE];
	[img saveAsKOImage:@"test_psf_quad_k.img"];

	[img setConst:1.0];
    [img sinFilt:1 ang: 0.0];
	[img saveAsKOImage:@"test_psf_quad_filter.img"];
}

void
make_ref(RecImage *f, float cy, float cx, float a, float th)
{
	int		i, j, ix;
	float	cs, sn, tth, r, x, y, re, im;
	int		n = [f xDim];
	float	*p = [f data];
	float	*q = p + [f dataLength];

	ix = (cy + n/2) * n + cx + n/2;
	p[ix] += a * cos(-th);
	q[ix] += a * sin(-th);
}

RecImage *
mkFilt(int n, BOOL cosFlg)
{
	RecImage	*f;
	int			i, j, ix;
	float		x, y, th;
	float		*p;

	f = [RecImage imageOfType:RECIMAGE_REAL xDim:n yDim:n];
	p = [f data];

	for (i = 0; i < n; i++) {
		y = (i - n/2);
		for (j = 0; j < n; j++) {
			x = (j - n/2);
            ix = i * n + j;
			th = atan2(y, x);
			if (cosFlg) {
				p[ix] = cos(th);
			} else {
				p[ix] = sin(th);
			}
		}
	}
	return f;
}

// SNR of end result
float
snr(RecImage *img)
{
    int     i, j, xDim;
    int     x, y;
    float   s, n;
    float   *p, *q;
    float   re, im, mg;

    p = [img real];
    q = [img imag];
    xDim = [img xDim];
    s = n = 0;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            x = j * 12 + 40;
            y = i * 12 + 40;
            re = p[y * xDim + x];
            im = q[y * xDim + x];
            mg = sqrt(re*re + im*im);
            s += mg;
        }
    }
    s /= (5 * 5);

    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            x = j;
            y = i;
            re = p[y * xDim + x];
            im = q[y * xDim + x];
            mg = sqrt(re*re + im*im);
            n += mg;
        }
    }
    n /= (32 * 32);

    return s / n;
}

void
dumpSig(RecImage *img1, RecImage *img2)
{
    FILE    *fp;
    int     i, j, xDim;
    int     ix;
    int     x, y;
    float   s, n;
    float   *p1, *q1;
    float   *p2, *q2;
    float   re, im, mg1, mg2;

    fp = fopen("PH_sig.txt", "w");
    p1 = [img1 real]; q1 = [img1 imag];
    p2 = [img2 real]; q2 = [img2 imag];
    xDim = [img1 xDim];
    s = n = 0;
    ix = 0;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++, ix++) {
            x = j * 12 + 40;
            y = i * 12 + 40;
            re = p1[y * xDim + x];
            im = q1[y * xDim + x];
            mg1 = sqrt(re*re + im*im);
            re = p2[y * xDim + x];
            im = q2[y * xDim + x];
            mg2 = sqrt(re*re + im*im);

            fprintf(fp, "%d %f %f\n", ix, mg1, mg2);
        }
    }
    fclose(fp);
}

void
calc_cnr_ft(RecImage *f)
{
	float	*p = [f data] + [f dataLength];	// imag part
	int		xdim = [f xDim];
	int		i, n = xdim * 32;
//	int		ix = xdim * 40 + 39;
	int		ix = xdim * 63 + 64;
	float	mn, sd, pk;

	printf("ft\n");

	Num_avevar(p, n, &mn, &sd);
	sd = sqrt(sd);
	pk = p[ix];
	printf("mean = %f, var = %f, pk = %f, cnr = %f\n", mn, sd, pk, pk/sd);
}

void
calc_cnr_kflt(RecImage *f)
{
	float	*p = [f data];	// real part
	int		xdim = [f xDim];
	int		i, n = xdim * 32;
//	int		ix = xdim * 40 + 40;
	int		ix = xdim * 64 + 64;
	float	mn, sd, pk;

	printf("kFlt\n");

	Num_avevar(p, n, &mn, &sd);
	sd = sqrt(sd);
	pk = p[ix];
	printf("mean = %f, var = %f, pk = %f, cnr = %f\n", mn, sd, pk, pk/sd);
}

void
calc_angle(RecImage *f)
{
	float	*p = [f data];	// real part
	float	*q = p + [f dataLength];
	int		xdim = [f xDim];
	int		i, j, k, n = 25;
	int		ix = xdim * 40 + 40;
	float	re, im, mg, phs, true_angle;

	printf("calc_angle\n");

	k = 0;
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++, k++) {
			true_angle = (k - n/2) * 2 * M_PI * rng / n;
			ix = i * 12 * xdim + j * 12 + xdim * 40 + 40;
			re = p[ix];
			im = q[ix];
			mg = re * re + im * im;
			mg = sqrt(mg);
			phs = -atan2(im, re) * M_PI * rng / M_PI;
			printf("%f %f\n", true_angle, phs);
		}
	}
}

//
void
kSpEn(RecImage *f)
{
    RecLoop     *kx, *ky;
    RecImage    *g;
    RecImage    *param;
    NSPoint     pos;

    kx = [RecLoop loopWithDataLength:256];
    ky = [RecLoop loopWithDataLength:256];
    [f replaceLoop:[f xLoop] withLoop:kx];
    [f replaceLoop:[f yLoop] withLoop:ky];
    [f conjugate];

    g = [f copy];
    [g magnitude];
    pos = [g findPeak2D];
    pos.x -= [g xDim]/2;
    pos.y -= [g yDim]/2;
    printf("x/y = %f/%f\n", pos.x, pos.y);

    param = [RecImage pointPoint:pos];
    f = [f ftShiftBy:param];
    [f saveAsKOImage:@"KSE_raw.img"];
// ft
    g = [f copy];
    [g fft2d:REC_FORWARD];
    [g saveAsKOImage:@"KSE.img"];
// p-sprctrum
    g = pSpectPt(f, 84, 141, YES);     // noise
    [g saveAsKOImage:@"KSE_raw.img"];
// p-space filter
    g = [f copy];
    g = pCirc(f, 32, YES);    // p-filt
    [g magnitude];
//    [g takeRealPart];
    [g saveAsKOImage:@"KSE_pf.img"];
    g = pFilt(g, 2);
    [g saveAsKOImage:@"KSE_pfi.img"];
}

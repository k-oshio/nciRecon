//
//	nciRec12 recon (mag)
//	7-7-2018		forked off from nciRec9
//	this is for NIRS experiments
//	new scheme using GR-EPI and long delay
//	second phase only (first phsae is done by nciRec11) ###

//	=== plans ===
//  1. make each proc function (keep main func short) ###
//	2. phasor analysis
//		calc phasor in one step, and save result		11-10-2018
//		process phasor for different position / block (in this program. fast enough)
//	3. BOLD again
//		* average over 1 sec -> 380 pts
//		average for one set (20 sec, 200 pts), and take PCA
//	== stop and phase again -> nciRec13
//
// ======= next phase (after QST final report) ========
//	1. LPF -> X
//  2. ICA for every mode (phs too)     7-12-2020


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

void		time_pca(RecImage *img);
void        time_ica(RecImage *img, int nComp);
RecImage    *time_phasor(RecImage *img);
RecImage    *calc_resp(RecImage *img, int nComp);
void        mov_avg(RecImage *img, int x, int y);
RecImage    *kern_avg(RecImage *img, int dly);
void		mk_kern(void);			// Gabor kernel
void		test_model(void);			// q & d model
RecImage	*phasor_filt(RecImage *img, float *re, float *im, int n);
RecImage	*phasor_hist(RecImage *img, int x, int y, float scale);
NSString	*base, *work, *path;
RecImage	*tmp_img;

int         dda = 100;    // 10sec * 10phs

// phasor
float		kc[11], ks[11];
float		frq = 0.8;
float		wd = 0.5;
int			nblk = 4;

int         study, series;
int         single_study;

// block average
int         bSize;
int			avgMode;	// 0: no avg (hires), 1: low res (1 fr/sec)
                        // 2: hires, 20sec win, 3: hires, 20sec win, sp mask
                        // 4: hires, 2sec win (for ser 1)
int         filter;     // 0: no flt, 1: gauss diff (freq), 2: sinusoid (time)

// ICA
int			nComp;

int
main(int ac, char *av[])
{
    RecImage        *img, *tmp_img;
	RecImage		*mgmask, *spmask, *avg;
	BOOL			bold;
    BOOL            mask_on;
	RecLoop			*tmLp, *phsLp, *avgLp, *blkLp, *xLp, *yLp;
	int				xDim = 64, yDim = 64;
	int				tDim, iDim, oDim, nPixel;
    int             i, j, k;
    float           *p, *pp;

	// PCA
//	RecImage		*img_a, *img_u;
//	Num_mat			*A;
//	Num_svd_result	*sres;
	// PCA filter
//	RecImage		*pca, *slice, *err;
//	float			*p, *q, *pp, *qq;
//	int				i, j, k, slc, blk;
	
    printf("nciRec12 (7-7-2018)\n");

    @autoreleasepool {
		img = nil;

// === test model
	if (0) {
		test_model();
		exit(0);
	}

// === read raw ===
//      3-2, *3-3, 4-2, *4-3, 6-2, *6-3
// ===== study loop ====
single_study = 3;
series = 2;
for (study = single_study; study <= single_study; study++) {
//for (study = 1; study <= 7; study++) {
    if (study == 5) continue; // skip
			//   1 : 2018-6-11-1
			//   2 : 2018-6-11-2 (V1 center)
			//  *3 : 2018-6-27
			//  *4 : 2018-6-29
			//   5 : 2019-1-21 (skip)
			//  *6 : 2019-2-4
			//   7 : 2019-2-15
			//==== 1.0 sec protocol
			// 1 : 10 * 10 dda, 10 * 380 (NCI, 1 stim/sec)
			// 2 : 10 * 10 dda, 10 * 380 (BOLD with short stim) 
			// 3 : 10 * 10 dda, 10 * 380 (BOLD with continuous stim)
			//==== 1.2 sec protocol
			// 2,3: 10 * 10 dda, 12 * 320 
		bold = YES;
                //         dwSize, blkSize
		avgMode = 5;	// 0: 1, 3200, 1: 10, 3200, 2: 1, 200, 240 3: 1, 200, sp mask 4: 1, 20
		switch (avgMode) {
        case 0 :
            nComp = 10;
            mask_on = YES;
            break;
        case 1 :
            nComp = 20;
            mask_on = YES;
            break;
		case 2 :
			nComp = 20;
			mask_on = NO;
			break;
        case 3:
            nComp = 20;
            mask_on = YES;
            break;
        case 4:
            nComp = 10;
            mask_on = YES;
            break;
        case 5:
            nComp = 20;
            mask_on = YES;
            break;
		}
  
        printf("stu:%d, ser:%d, avg:%d\n", study, series, avgMode);

	printf("read raw\n");
		switch (study) {
		case 1:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-1";
			break;	// study 1
		case 2:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-2";
			break;	// study 2
		case 3:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0627";
			break;	// study 3
		case 4:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0629";
			break;	// study 4
		case 5:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0121";
			break;	// study 5
		case 6:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0204";
			break;	// study 6
		case 7:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0215";
			break;	// study 7
		}
        work = [NSString stringWithFormat:@"%@/results/%d", base, series];
		path = [NSString stringWithFormat:@"%@/IMG_comb", work];
		img = [RecImage imageWithKOImage:path]; // input image

		avg = [img avgForLoop:[img zLoop]];
		path = [NSString stringWithFormat:@"%@/IMG_avg", work];
		[avg saveAsKOImage:path];

	// sp mask
		spmask = [avg copy];
		[spmask magnitude];
		p = [spmask data];
		iDim = [spmask dataLength];
		for (i = 0; i < iDim; i++) {
			if (i < iDim*32/64) {                      // 32
				p[i] = 0;
			} else {
				p[i] = 1;
			}
		}
	// mag mask
		mgmask = [avg copy];
		[mgmask magnitude];
		[mgmask thresAt:0.2];
        if (mask_on) {
            [mgmask multByImage:spmask];
        }
	//	[mgmask hGauss2DLP:0.5];	// "Half" gaussian filter
		p = [mgmask data];
		nPixel = 0;
		for (i = 0; i < [mgmask dataLength]; i++) {
			if (p[i] > 0) {
				nPixel++;
			}
		}
//		printf("nPixel = %d\n", nPixel);

// ===== block ===== (current 11-10-2018)
//		if (series == 2 || series == 3) {
        // select sub-set for processing
        xLp = [RecLoop loopWithDataLength:xDim];
        yLp = [RecLoop loopWithDataLength:yDim];
        tDim = 3600; //1200;	// 3800 total, (100-on / 100-off) x 16
        if (study > 5) {
            tDim = 3360;
        }
        iDim = xDim * yDim;
    //	avgLp = [img crop:avgLp to:380 startAt:dda];
   //     tmLp = [img crop:[img zLoop] to:tDim startAt:dda + 0];
        switch (avgMode) {
        case 0 :
        case 1 :
            tDim += dda;
            tmLp = [img crop:[img zLoop] to:tDim startAt:0];
            break;
        case 2:
        case 3:
            tmLp = [img crop:[img zLoop] to:tDim startAt:dda];
            break;
        case 4:
            tDim = 1800;
            tmLp = [img crop:[img zLoop] to:tDim startAt:dda + 600];
            break;
        case 5 :
            tDim += dda;
            tmLp = [img crop:[img zLoop] to:tDim startAt:0];
            break;
        }

        [img multByImage:mgmask];

        // copied from nciRec13
        if (bold) {
            [img magnitude];
         //       [img baseline];
            path = [NSString stringWithFormat:@"%@/IMG_mag", work];
            [img saveAsKOImage:path];
        } else {
            avg = [img avgForLoop:[img zLoop]];
            [avg phase];
            [img phase];        // need to remove baseline phase first ...
        //    avg = [img avgForLoop:[img zLoop]];
            [img subImage:avg];
            [img unwrap0d];
            path = [NSString stringWithFormat:@"%@/IMG_phs", work];
            [img saveAsKOImage:path];
        }

    // baseline
        avg = [img avgForLoop:[img zLoop]]; // mean of all
        
        [img subImage:avg];
        path = [NSString stringWithFormat:@"%@/IMG_img0", work];
        [img saveAsKOImage:path];
        
    // filter
        switch (filter) {
        case 0:
            break;
        case 1:
            [img gaussDiff:0.5 forLoop:[img zLoop]];
            break;
        case 2:
            [img fSinWLen:10 phase:0.5 forLoop:[img zLoop]];
            break;
        }
path = [NSString stringWithFormat:@"%@/IMG_tflt", work];
[img saveAsKOImage:path];
exit(0);


        switch (avgMode) {
        case 0 :
        // no average
            bSize = 1;
            tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
            [tmp_img copyImageData:img];
            img = [tmp_img avgForLoop:[tmp_img zLoop]];
            break;
        case 1 :
        // time average (1 frame/sec)
            bSize = 10;
            tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
            [tmp_img copyImageData:img];
            img = [tmp_img avgForLoop:[tmp_img zLoop]];
            break;
        case 2 :    // sp mask off
        case 3 :    // sp mask on
           // block avg (10 frames/sec, for 20sec)
            if (study <= 5) {
                bSize = 200;
            } else {
                bSize = 240;
            }
            tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
            [tmp_img copyImageData:img];
            img = [tmp_img avgForLoop:[tmp_img topLoop]];
            break;
        case 4:     // 2sec win
            if (study <= 5) {
                bSize = 40;
            } else {
                bSize = 48;
            }
            tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
            [tmp_img copyImageData:img];
        //    mov_avg(img, 26, 43);   // off
            mov_avg(img, 23, 43); // on
            img = [tmp_img avgForLoop:[tmp_img topLoop]];
            break;
        case 5 :
        // time average (1 frame/sec), with response kernel
            bSize = 10;
            tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
            [tmp_img copyImageData:img];
        //    img = [tmp_img avgForLoop:[tmp_img zLoop]];
            img = kern_avg(tmp_img, 1); // input, delay (samples)
            [img dumpLoops];
            break;
        }
        if (bold) {
            work = [NSString stringWithFormat:@"%@/results/%d/mag/%d", base, series, avgMode];
        } else {
            work = [NSString stringWithFormat:@"%@/results/%d/phs/%d", base, series, avgMode];
        }
//[img magnitude];
//[img phase];
        path = [NSString stringWithFormat:@"rm %@/IMG_*", work];
        system([path UTF8String]);
        path = [NSString stringWithFormat:@"%@/IMG_subsample", work];
        [img saveAsKOImage:path];

//		time_pca(img);
        time_ica(img, nComp);
//        tmp_img = time_phasor(img);
//        path = [NSString stringWithFormat:@"%@/IMG_psr", work];
//        [tmp_img saveAsKOImage:path];

            
    printf("study #%d done\n", study);
}   // ===== study loop ==========
    
    }   // @autoreleasepool

    return 0;
}

// add mask, path modifier
void
time_pca(RecImage *img)
{
	int				tDim, iDim;
	int				xDim, yDim;
	int				i;
	RecImage		*img_a, *img_u;
	RecImage		*pca;
	Num_mat			*A;
	Num_svd_result	*sres;

	tDim = [img zDim];
	xDim = [img xDim];
	yDim = [img yDim];
	iDim = xDim * yDim;
	img_a = [RecImage imageOfType:[img type] xDim:iDim yDim:tDim];
	pca = [RecImage imageOfType:[img type] xDim:xDim yDim:yDim zDim:tDim];
	
	[img_a copyImageData:img];
	[img_a dumpLoops];
	A = Num_im_to_m(img_a);
//path = [NSString stringWithFormat:@"%@/IMG_A", work];
//[img_a saveAsKOImage:path];
//saveAsKOImage(A, path);
	sres = Num_svd(A);
	img_a = Num_m_to_im(sres->Vt);
	[pca copyImageData:img_a];
	path = [NSString stringWithFormat:@"%@/IMG_PCA", work];
	[pca saveAsKOImage:path];
	img_u = Num_m_to_im(sres->U);
	[img_u trans];
	[img_u crop:[img_u yLoop] to:32 startAt:0];
	path = [NSString stringWithFormat:@"%@/IMG_U", work];
	[img_u saveAsKOImage:path];
	for (i = 0; i < 20; i++) {
		printf("%d %f\n", i, sres->s->data[i]);
	}
	tmp_img = [img_u copy];
	[tmp_img fft1d:[img_u xLoop] direction:REC_INVERSE];
	path = [NSString stringWithFormat:@"%@/IMG_U_ft", work];
	[tmp_img saveAsKOImage:path];
}

NumMatrix *
makeRespMat(int nt) // time dim
{
    NumMatrix   *mat;
    float       *p;
    float       hi, lo;
    int         i, j, k, del, dstp; // 1 sec interval
    int         nd;
    int         nstim = 1;
    int         non, noff, blklen;
    // global : study, series, avgMode

//    if (study <= 5) {
//        nd = 20;
//    } else {
//        nd = 24;
//    }
    nd = 1;
    switch (avgMode) {
    case 1 :
        mat = [NumMatrix matrixOfType:NUM_REAL nRow:nt nCol:nd];
        dstp = 1;

        if (study <= 5) { // pass stmlen later
            blklen = 20;
        } else {
            blklen = 24;
        }
        non = blklen/2;
        noff = blklen - non;
        nstim = nt / blklen;
        hi = 1.0;
//        lo = -1.0 / noff;
        lo = 0;
        for (k = 0; k < nd; k++) {
            del = dstp * k;
            p = [mat data] + nt * k;
            for (j = 0; j < dda/10; j++) {
                p[j] = lo;
            }
            for (i = 0; i < nstim; i++) {
                p = [mat data] + nt * k + i * blklen + dda/10;
                for (j = del; j < non + del; j++) {
                    p[j] = (float)(j - del) * hi / non + lo;
                }
                for (; j < blklen + del; j++) {
                    p[j] = (float)(blklen - j - del) * hi / noff + lo;
                }
            }
        }
        break;
    case 2:
        mat = [NumMatrix matrixOfType:NUM_REAL nRow:nt nCol:nd];
        blklen = nt;
        if (series == 1) {
            dstp = 1;
            non = blklen/2;
            noff = blklen - non;
        } else {
            dstp = 10;
            non = blklen/2;
        }
        nstim = nt / blklen;
        for (k = 0; k < nd; k++) {
            del = dstp * k;
            p = [mat data] + nt * k;
            for (j = 0; j < nt; j++) {
                p[j] = -0.5;
            }
            for (i = 0; i < nstim; i++) {
                p = [mat data] + nt * k + i * blklen;
                for (j = del; j < blklen + del; j++) {
                    p[j] = 0.5;
                }
            }
        }
        break;
    }

    return mat;
}

// ### started on 12-27-2018 (copied pca)
// ### 7-12 probably ok now...
void
time_ica(RecImage *img, int nComp)
{
    int             tDim, iDim;
    int             xDim, yDim;
//    int             i, j;
    RecImage        *img_a, *img_sp, *img_tm;
    RecImage        *pca, *pcat;
    NumMatrix       *A, *tmp, *sg;
    NSDictionary    *res;
    RecImage        *bold, *nci, *at, *col;
    NSString        *path;
    NumMatrix       *Resp, *P, *S, *C;
    float           *buf, *p, *s, sum;       // -> make this Matrix obj
    // global : study, series, avgMode

path = [NSString stringWithFormat:@"%@/IMG_in", work];
[img saveAsKOImage:path];

    xDim = [img xDim];
    yDim = [img yDim];
    tDim = [img zDim];
    iDim = xDim * yDim;
    img_a = [RecImage imageOfType:RECIMAGE_REAL xDim:iDim yDim:tDim];   // this is correct (trans of doc)
    [img_a copyImageData:img];
    A = [img_a toMatrix];
path = [NSString stringWithFormat:@"%@/IMG_A", work];
[A saveAsKOImage:path];

    res = [A icaForNC:nComp];
    img_sp = [img copy];
    [img_sp crop:[img_sp zLoop] to:nComp startAt:0];
// *PCA (U)
    [img_sp copyImageData:[[res objectForKey:@"U"] toRecImage]];
    path = [NSString stringWithFormat:@"%@/IMG_PCA", work];
    [img_sp saveAsKOImage:path];
// *PCAt (Vt)
    img_tm = [[res objectForKey:@"Vt"] toRecImage];
    path = [NSString stringWithFormat:@"%@/IMG_PCAt", work];
    [img_tm saveAsKOImage:path];
// PCA sg
    sg = [res objectForKey:@"Sg"];
    path = [NSString stringWithFormat:@"%@/IMG_Sg", work];
    [tmp saveAsKOImage:path];
// *ICA (Y) 
    [img_sp copyImageData:[[res objectForKey:@"Y"] toRecImage]]; // -> make "copyMatrixData" method
    path = [NSString stringWithFormat:@"%@/IMG_ICA", work];
    [img_sp saveAsKOImage:path];
// *color map
    bold = [img_sp sliceAtIndex:0];
    nci  = [img_sp sliceAtIndex:3];
    at   = [img_sp sliceAtIndex:2];
    [bold magnitude];
    [nci magnitude];
    [at magnitude];
    col = [img_sp makeColorWithR:nci G:bold B:at];
    [col saveAsKOImage:@"IMG_col"];
// *ICAt (WX)
    img_tm = [[res objectForKey:@"WX"] toRecImage];
    [img_tm trans];
    path = [NSString stringWithFormat:@"%@/IMG_ICAt", work];
    [img_tm saveAsKOImage:path];    // WX
// *ICA W
    img_tm = [[res objectForKey:@"W"] toRecImage];
    path = [NSString stringWithFormat:@"%@/IMG_W", work];
    [img_tm saveAsKOImage:path];    // W
// *ICA sg
    tmp = [res objectForKey:@"WSg"];
    path = [NSString stringWithFormat:@"%@/IMG_WSg", work];
    [tmp saveAsKOImage:path];

// ===== projection onto stim vector -> X
    Resp = makeRespMat(tDim); // time dim, delay dim
    path = [NSString stringWithFormat:@"%@/IMG_Resp", work];
    [Resp saveAsKOImage:path];

    P = [res objectForKey:@"Vt"];
    S = [[res objectForKey:@"Sg"] diagMat];
    C = [[Resp trans] multByMat:[P trans]];
    tmp = [C multByMat:P];
    path = [NSString stringWithFormat:@"%@/IMG_Qt", work];
    [tmp saveAsKOImage:path];

    [img_sp crop:[img_sp zLoop] to:[Resp nCol] startAt:0];
    tmp = [res objectForKey:@"U"];
    tmp = [S multByMat:tmp];
    tmp = [C multByMat:tmp];
    path = [NSString stringWithFormat:@"%@/IMG_Qsp", work];
    [img_sp copyImageData:[tmp toRecImage]];
    [img_sp saveAsKOImage:path];

[img_sp clear];
    tmp = [[Resp trans] multByMat:A]; // ### dim not correct
    [img_sp copyImageData:[tmp toRecImage]];
//    [img_sp subImage:[img_sp sliceAtIndex:0]];
    path = [NSString stringWithFormat:@"%@/IMG_BOLD", work];
    [img_sp saveAsKOImage:path];
    


}

// visualize delay
RecImage *
time_phasor(RecImage *img)
{
    int         i, j;
    float       *pp;
    float       *p, *q;
    int         tDim, iDim;
    float       wLen;
    float       *cs, *sn, th, sig;
    float       sumr, sumi;
    RecImage    *psr;

    tDim = [img zDim];
    iDim = [img xDim] * [img yDim];

    psr = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:[img xDim] yDim:[img yDim]];
    switch (avgMode) {
    case 1 :
        if (study <= 5) { // pass stmlen later
            wLen = 20;
        } else {
            wLen = 24;
        }
        break;
    case 2 :
        wLen = tDim;
        break;
    }
printf("wlen = %f\n", wLen);

    cs = (float *)malloc(sizeof(float) * tDim);
    sn = (float *)malloc(sizeof(float) * tDim);
    for (i = 0; i < tDim; i++) {
        th = (float)i * M_PI * 2 / wLen;
        cs[i] = cos(th);
        sn[i] = sin(th);
        printf("%d %f %f\n", i, cs[i], sn[i]);
    }
    
    pp = [img data];
    p = [psr real];
    q = [psr imag];
    for (i = 0; i < iDim; i++) {
        sumr = sumi = 0;
        for (j = 0; j < tDim; j++) {
            sig = pp[j * iDim + i];
            sumr += sig * cs[j];
            sumi += sig * sn[j];
        }
        p[i] = sumr;
        q[i] = sumi;
    }

    return psr;
}

// phasor filter
RecImage *
phasor_filt(RecImage *img, float *re, float *im, int nk)
{
	RecImage	*flt;
	float		*pp, *p, *q;
	float		sumr, sumi, th, cs, sn;
	int			i, j, k, jj, nimg, npix;

	nimg = [img zDim];
	flt = [RecImage imageOfType:RECIMAGE_COMPLEX withImage:img];
	p = [flt real];
	q = [flt imag];
	pp = [img data];
	npix = [img xDim] * [img yDim];

	for (i = 0; i < npix; i++) {
//	printf("%d %f\n", i, th);
		for (j = 0; j < nimg; j++) {
			th = (float)j * frq * M_PI;
			cs = cos(-th);
			sn = sin(-th);
			sumr = sumi = 0;
			for (k = 0; k < nk; k++) {
				jj = j + k - nk/2;
				if (jj < 0) continue;
				if (jj >= nimg) continue;
				sumr += pp[jj * npix + i] * re[k];
				sumi += pp[jj * npix + i] * im[k];
			}
			// ### phase correction
			p[j * npix + i] = cs * sumr + sn * sumi;
			q[j * npix + i] = sn * sumr - cs * sumi;
		}
	}
	return flt;
}

// ### include block avg code ###
// --> separate input into phases beforehand img:[phs, avg, yLp, xLp]
RecImage *
phasor_hist(RecImage *img, int posx, int posy, float scale)
{
	RecImage	*hst;
	float		*p, *q;
	float		*pp;
	int			i, k, len = [img zDim];
	int			xDim = [img xDim];
	int			yDim = [img yDim];
	int			nBlk = [img topDim];
	float		x, y;
	int			ix, iy;
	int			dim = 256;

	hst = [RecImage imageOfType:RECIMAGE_REAL xDim:dim yDim:dim zDim:nBlk];

//printf("nBlk = %d, len = %d, xDim = %d\n", nBlk, len, xDim);
	for (k = 0; k < nBlk; k++) {
		p = [img real] + k * len * xDim * yDim;
		q = [img imag] + k * len * xDim * yDim;
		pp = [hst data] + k * dim * dim;
		for (i = 0; i < len; i++) {
			x = p[i * xDim * yDim + posy * xDim + posx];
			y = q[i * xDim * yDim + posy * xDim + posx];
			ix =  (int)(x * scale) + dim/2;
			iy = -(int)(y * scale) + dim/2;
			if (ix < 0) ix = 0;
			if (ix >= dim) ix = dim - 1;
			if (iy < 0) iy = 0;
			if (iy >= dim) iy = dim - 1;
			pp[iy * dim + ix] += 1;
		}
	}
	pp = [hst data];
	for (i = 0; i < [hst dataLength]; i++) {
		pp[i] = log(pp[i] + 1) * 100;
	}

	return hst;
}

void
mk_kern(void)
{
	int		i, n = 11;
	float	x, w, sum, re, im;
//	float	c[50], s[50];

	sum = re = im = 0;
	for (i = 0; i < n; i++) {
		x = ((float)i - (float)n/2) * frq * M_PI;
	printf("%d %f\n", i, x);
		w = exp(-x*x/wd/wd/n/n);	// make this spline
	//	kc[i] = cos(x) * w;
	//	ks[i] = sin(x) * w;
		kc[i] = cos(x + M_PI/4) * w;
		ks[i] = cos(x - M_PI/4) * w;
		sum += kc[i] * ks[i];
		re += kc[i] * kc[i];
		im += ks[i] * ks[i];
	//	printf("%d %f %f %f\n", i, c[i], s[i], w);
	}
	re = sqrt(re);
	im = sqrt(im);

	for (i = 0; i < n; i++) {
		kc[i] /= re;
		ks[i] /= im;
//		printf("%d %f %f\n", i, c[i], s[i]);
	}
	printf("proj = %f %f %f\n", sum, re, im);

	re = im = 0;
	for (i = 0; i < n; i++) {
		re += kc[i];
		im += ks[i];
	}
	printf("re = %f, im = %f\n", re, im);
//	kc[n/2] -= re;
	re = im = 0;
	for (i = 0; i < n; i++) {
		re += kc[i];
		im += ks[i];
		printf("%d %f %f\n", i, kc[i], ks[i]);
	}
	printf("re = %f, im = %f\n", re, im);


}

void
test_model(void)			// q & d model
{
	int		i;
	float	x, a;
	float	pd = 6500;
	float	d1 = 2.5;
	float	d2 = 3.5;
	float	t1 = 5.0;
	float	t2 = 2.0;

	for (i = 0; i < 200; i++) {	// x 10msec
		x = (float)i / 10;	// sec
		if (x >= 0 && x < d1) {
			a = 0;
		} else
		if (x < 10 + d2) {
			a = 1.0 - exp(-(x - d1) / t1);
		} else {
			a = 1.0 - exp(-(10.0 + d2 - d1) / t1);
			a *= exp(-(x - 10 - d2) / t2);
		}
		a = pd * (a - 0.4);
		printf("%d %f\n", i, a);
	}
}

RecImage *
calc_resp(RecImage *img, int nComp)
{
    RecImage    *ir, *stim, *stmg, *stphs, *tmg, *tphs;
    RecLoopControl  *lc;
    int         i, j, ix, n, len, skip;
    float       *p1, *p2, *q1, *q2;
    float       nz = 0.01;

    ir = [img copy];
    [ir fft1d:[ir zLoop] direction:REC_INVERSE];
    path = [NSString stringWithFormat:@"%@/IMG_ft", work];
    [ir saveAsKOImage:path];
    tmg = [ir copy];
    [tmg magnitude];
    tphs = [ir copy];
    [tphs phase];

    stim = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:[img zLoop], nil];
    p1 = [stim data];
    n = [img zDim];
    len = [img xDim] * [img yDim];
    skip = [img skipSizeForLoop:[img zLoop]];
    for (j = 0; j < n/2; j++) {
        p1[j] = 0.5;
        p1[j + n/2] = -0.5;
    }
    [stim fft1d:[stim xLoop] direction:REC_INVERSE];
//    path = [NSString stringWithFormat:@"%@/IMG_stim", work];
//    [stim saveAsKOImage:path];
    stmg = [stim copy];
    [stmg magnitude];
    [stmg scaleToVal:1.0];
    stphs = [stim copy];
    [stphs phase];
    p2 = [stmg data];
    q2 = [stphs data];

//    [tmp_img cpxDivImage:stim];    // doesn't work -> do manually
    lc = [img control];
    [lc deactivateLoop:[img zLoop]];
    len = [lc loopLength];
    for (i = 0; i < len; i++) {
        p1 = [tmg currentDataWithControl:lc];
        q1 = [tphs currentDataWithControl:lc];
        for (j = ix = 0; j < n; j++, ix+=skip) {
            if (p2[j] < nz || abs(j - n/2) > nComp) {
                p1[ix] = 0;
                q1[ix] = 0;
            } else {
                p1[ix] /= p2[j] + nz;
                q1[ix] += q2[j];
            }
        }
        [lc increment];
    }

    ir = [tmg copy];
    [ir makeComplexWithPhs:tphs];

    [ir fft1d:[ir zLoop] direction:REC_FORWARD];
//    path = [NSString stringWithFormat:@"%@/IMG_ift", work];
//    [ir saveAsKOImage:path];
    [ir takeRealPart];

    return ir;
}

void
mov_avg(RecImage *img, int x, int y)
{
    int         i, j, n, kwin;
    int         xDim, yDim, imgSize, ofs;
    float       *bufi, *p;
    float       *bufre, *bufim;
    RecImage    *res;
    float       *cs, *sn, th, re, im;
    int         ncyc = 1;

//    [img dumpLoops]; // 90, 40, 64, 64

    xDim = [img xDim];
    yDim = [img yDim];
    imgSize = xDim * yDim;
    kwin = 10 * ncyc;
    cs = (float *)malloc(sizeof(float) * kwin);
    sn = (float *)malloc(sizeof(float) * kwin);
    for (i = 0; i < kwin; i++) {
        th = (float)i * ncyc * 2 * M_PI / kwin;
        cs[i] = cos(th);
        sn[i] = sin(th);
    //    printf("%d %f %f\n", i, cs[i], sn[i]);
    }
    n = [img zDim];

    res = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n/kwin];
    bufre = [res real];
    bufim = [res imag];
    bufi = (float *)malloc(sizeof(float) * n * kwin);

    p = [img data];
    ofs = y * xDim + x;
    for (i = 0; i < n; i++) {
        re = 0;
        re += p[i * imgSize + ofs];
        re += p[i * imgSize + ofs + 1];
        re += p[i * imgSize + ofs + xDim];
        re += p[i * imgSize + ofs + xDim + 1];
        bufi[i] = re / 4;
    //    printf("%d %f\n", i, re);
    }
    
    for (i = 0; i < n/kwin; i++) {
        re = im = 0;
        for (j = 0; j < kwin; j++) {
            re += cs[j] * bufi[i * kwin + j];
            im += sn[j] * bufi[i * kwin + j];
        }
        bufre[i] = re / kwin;
        bufim[i] = im / kwin;
        //    printf("%d %f %f\n", i, bufm[i], bufp[i]);
         //   printf("%d %f %f\n", i, re, im);
    }
    [res gauss1DLP:0.02 forLoop:[res xLoop]];
    [res saveAsKOImage:@"IMG_res"];
        


}

RecImage *
kern_avg(RecImage *img, int dly)
{
    RecImage        *res;
    RecLoop         *lp;
    RecLoopControl  *lc;
    float           *p, *q, val;
    int             i, j;
    int             n, skip, len, imgSize;
    
    lc = [img control];
    lp = [img zLoop];
    [lc removeLoop:lp];
    res = [RecImage imageOfType:RECIMAGE_REAL withControl:lc];
    skip = [img skipSizeForLoop:lp];
    len = [lp dataLength];
    imgSize = [img xDim] * [img yDim];
// mult by kern
    [lc rewind];
    n = [lc loopLength];
printf("n = %d, skip = %d, len = %d imgSize = %d\n", n, skip, len, imgSize);
    q = [res data];
    for (i = 0; i < n; i++) {
        p = [img currentDataWithControl:lc];    // each pixel
        val = 0;
        for (j = 0; j < len; j++) {
            if (j < dly || j >= dly + len/2) {
                val -= p[j*skip];
            } else {
                val += p[j*skip];
            }
        //    val += p[j*skip];
        }
        q[i] = val / len;
        [lc increment];
    }
//res = [img avgForLoop:lp];
    return res;
}



//
//	nciRec13 recon (phase)
//	second phase only (first phsae is done by nciRec11) ###

//	=== plans ===
//	4. phase again, forked off from nciRec12 (BOLD)
//		*take dif along time
//		*blk average
//		1 Hz filter -> 8 Hz is aliased to 2 Hz, with 10 Hz sampling -> don't remove 2 Hz
//		pca filter
//		correlation with input function (8 Hz flash, 2 cycles)
//	5. make directory structure to reserve results (6 toal)
//		mag / phase
//		low res / high res / freq domain (for cardiac component)
//	6. chk 0.83　Hz response with 1.0 Hz filter
//
// ======= next phase (after QST final report) ========
//	1. LPF / accum filter
//		remove 0th and 1st from delay loop
//		avgMode 3: ## 4: ##
	


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

void		time_pca(RecImage *img, int nflt);
void		time_ica(RecImage *img, int nComp);
void		mk_kern(void);			// Gabor kernel
RecImage	*phasor_filt(RecImage *img, float *re, float *im, int n);
RecImage	*phasor_hist(RecImage *img, int x, int y, float scale);
NSString	*base, *work, *path;
RecImage	*tmp_img;

// phasor
float		kc[11], ks[11];
float		frq = 0.8;
float		wd = 0.5;
int			nblk = 4;

// block average
int			avgMode;	// 0: low-res, no block, 1: hi-res, 20 sec block, 2: hi-res, 2 sec block, 3: on-off -> freq analysis
int			bSize;

// ICA
int			nComp;

int
main(int ac, char *av[])
{
    RecImage        *img, *bimg, *aimg;	// input, block, block average
	RecImage		*tmp_img;
	RecImage		*resp, *rindex;
	RecImage		*mgmask, *spmask, *roi, *avg, *var;
	int				study, series;
	RecLoop			*tmLp, *xLp, *yLp, *bLp, *avgLp;
	int				dda = 100;	// 10sec * 10phs
	int				xDim = 64, yDim = 64;
	int				tDim, iDim, nPixel;
	BOOL			bold;
	BOOL			mask_on;
	BOOL			dx;
	BOOL			perFilt;
	// PCA filter
//	RecImage		*pca, *slice, *err;
	float			*p, *pp;
	int				i, j, k;
	
    printf("nciRec13 (3-20-2019)\n");

    @autoreleasepool {
		img = nil;
	//	system("rm IMG_*");

// === tmp test
if (0) {
	RecImage	*img = [RecImage imageOfType:RECIMAGE_REAL xDim:64];
	float		*p = [img data];

	p[0] = 1.0;
	[img leakyIntWithTau:10.0 cutOff:0.05 forLoop:[img xLoop]];
	[img saveAsKOImage:@"integ.img"];

	exit(0);
}

// === test model
		if (0) {
			RecImage	*resp;
			RecLoop		*lp = [RecLoop loopWithDataLength:200];
			resp = genResp(lp, 0.1);
			exit(0);
		}
// === freq of stim signal
		if (0) {
			int			i, n = 240;
			RecImage	*st = [RecImage imageOfType:RECIMAGE_REAL xDim:n];
			float		*p = [st data];
			
			p[0] = 1;
			p[1] = -1;
			p[2] = 1;
			p[3] = -1;
			
			[st fft1d:[st xLoop] direction:REC_INVERSE];
			[st saveAsKOImage:@"stim.img"];
			exit(0);
		}
// === read raw ===
		study = 3;	// *: clean BOLD response
			//  1 : 2018-6-11-1		TR = 1.0s
			//  2 : 2018-6-11-2		TR = 1.0s
			// *3 : 2018-6-27		TR = 1.0s
			// *4 : 2018-6-29		TR = 1.0s
			// X5 : 2019-1-21		TR = 1.0s
			// *6 : 2019-2-4		TR = 1.2s
			//  7 : 2019-2-15		TR = 1.2s
		series = 2;
			//==== 1.0 sec protocol
			// 1 : 10 * 10 dda, 10 * 380 (NCI)
			// 2 : 10 * 10 dda, 10 * 380 (BOLD with short stim) 
			// 3 : 10 * 10 dda, 10 * 380 (BOLD with continuous stim)
			//==== 1.2 sec protocol
			// 2,3: 10 * 10 dda, 12 * 320 
		bold = YES;
		mask_on = YES;
		dx = NO;
		if (bold) {
			perFilt = NO;	// 1.0 Hz filter
		} else {
			perFilt = YES;	// 1.0 Hz filter
		}
		avgMode = 2;	// 0: full, 1: no block, hi/lo-res, 2: 20 or 24 sec (BOLD) block, *hi-res, *lo-res, response
						// 3: full, LPF at 4Hz, 6: on/off phasor
		nComp = 20;		// number of ICA components

	printf("read raw\n");
		switch (study) {
		case 1:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-1";
			break;
		case 2:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-2";
			break;
		case 3:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0627";
			break;
		case 4:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0629";
			break;
		case 5:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0121";
			break;
		case 6:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0204";
			break;
		case 7:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0215";
			break;
		}
	// input image
		path = [NSString stringWithFormat:@"%@/results/%d/IMG_comb", base, series];
		img = [RecImage imageWithKOImage:path];

	// output directories

		if (bold) {
			work = [NSString stringWithFormat:@"%@/results/%d/mag/%d", base, series, avgMode];
		} else {
			work = [NSString stringWithFormat:@"%@/results/%d/phs/%d", base, series, avgMode];
		}


	// initially no loop structure
		
	//	[img copyImageData:tmp_img];
		if (study == 3 || study == 4) {
			[img rotate:1];
		}

		avg = [img avgForLoop:[img zLoop]];
		path = [NSString stringWithFormat:@"%@/IMG_avg", work];
		[avg saveAsKOImage:path];

	// sp mask
		spmask = [avg copy];
		[spmask magnitude];
		p = [spmask data];
		iDim = [spmask dataLength];
		for (i = 0; i < iDim; i++) {
		//	if (i < iDim*36/64 || i > iDim*51/64) {
			if (i < iDim*35/64) {
				p[i] = 0;
			} else {
				p[i] = 1;
			}
		}
	// mag mask
		mgmask = [avg copy];
		[mgmask magnitude];
		[mgmask thresAt:0.25];
		if (mask_on) {
			[mgmask multByImage:spmask];		// sp mask on ######
		}
	//	[mgmask hGauss2DLP:0.5];	// "Half" gaussian filter
		p = [mgmask data];
		nPixel = 0;
		for (i = 0; i < [mgmask dataLength]; i++) {
			if (p[i] > 0) {
				nPixel++;
			}
		}
		printf("nPixel = %d\n", nPixel);

// ===== block ===== (current 11-10-2018)
			// select sub-set for processing
		//	xLp = [RecLoop loopWithDataLength:xDim];
		//	yLp = [RecLoop loopWithDataLength:yDim];
			xLp = [img xLoop];
			yLp = [img yLoop];
			tDim = 3600; //1200;	// 3800 total, (100-on / 100-off) x 18
			iDim = xDim * yDim;
		//	avgLp = [img crop:avgLp to:380 startAt:dda];
			tmLp = [img crop:[img zLoop] to:tDim startAt:dda + 0];

			[img multByImage:mgmask];

			if (bold) {
				[img magnitude];
				path = [NSString stringWithFormat:@"%@/IMG_mag", work];
				[img saveAsKOImage:path];
			} else {
				avg = [img avgForLoop:[img zLoop]];
				[avg phase];
				[img phase];		// need to remove baseline phase first ...
			//	avg = [img avgForLoop:[img zLoop]];
				[img subImage:avg];
				[img unwrap0d];
				path = [NSString stringWithFormat:@"%@/IMG_phs", work];
				[img saveAsKOImage:path];
	//			[img saveAsKOImage:path];
				
				if (dx) {
					img = [img dxForLoop:tmLp];		// differentiate along time ####
					path = [NSString stringWithFormat:@"%@/IMG_phs_d", work];
					[img saveAsKOImage:path];
				}
			}

			if (perFilt) {
		//		[img singleFreqFilt:[img zDim] * 0.2 forLoop:[img zLoop]];	// remove single frequency (single fourier coeff)
				[img cycFilt:10 forLoop:[img zLoop]];	// remove single frequency (single fourier coeff)
			}

			// ### calc bimg here (introduce loop structure into img)
			// block avg (10 frames/sec, for 20sec)
			if (study <= 5) {
				bSize = 200;	// 800
			} else {
				bSize = 240;
			}
			bLp = [RecLoop loopWithDataLength:bSize];
			avgLp = [RecLoop loopWithDataLength:[img zDim] / bSize];
			bimg = [RecImage imageOfType:[img type] withLoops:avgLp, bLp, yLp, xLp, nil];
			[bimg copyImageData:img];

			switch (avgMode) {
			case 0 :	// no proc (full PCA)
				time_pca(img, 2);
				break;
			case 1 :
			// time average (1 frame/sec)
			//	bSize = 20;
			//	bLp = [RecLoop loopWithDataLength:bSize];
			//	avgLp = [RecLoop loopWithDataLength:[img zDim] / bSize];
			//	bimg = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
			//	bimg = [RecImage imageOfType:[img type] withLoops:avgLp, bLp, yLp, xLp, nil];
			//	[bimg copyImageData:img];
				bimg = [bimg avgForLoop:avgLp];
				tmp_img = [bimg avgForLoop:[bimg zLoop]];
				[bimg subImage:tmp_img];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_avg", work];
				[bimg saveAsKOImage:path];
				break;
			case 2 :
				tmp_img = [bimg avgForLoop:avgLp];
				[tmp_img subImage:[tmp_img avgForLoop:bLp]];

				// BOLD (ok)
				[tmp_img gauss1DLP:0.08 forLoop:bLp];	// 0.08 BOLD
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_avg", work];
				[tmp_img saveAsKOImage:path];

				// projection to response function (BOLD & NCI) (spatial pattern of response)
				resp = genResp(bLp, 0.1);
				[tmp_img multBy1dImage:resp];
				tmp_img = [tmp_img avgForLoop:[tmp_img zLoop]];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_resp", work];
				[tmp_img saveAsKOImage:path];

				// ### select pixels with strongest response, and make 1D response func
				roi = [tmp_img copy];
				[roi multByImage:spmask];
				[roi thresAt:0.6];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_roi", work];
				[roi saveAsKOImage:path];
				
				tmp_img = [bimg copy];
				[tmp_img multByImage:roi];
				tmp_img = [tmp_img avgForLoop:xLp];
				tmp_img = [tmp_img avgForLoop:yLp];
				[tmp_img addConst:-[tmp_img meanVal]];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_respsig", work];
				[tmp_img saveAsKOImage:path];

				rindex = [tmp_img copy];

				// works... fix multBy1dImage:
				p = [rindex data];
				pp = [resp data];
				for (i = 0; i < [avgLp dataLength]; i++) {
					for (j = 0; j < bSize; j++) {
						p[i * bSize + j] *= pp[j];
					}
				}

			//	path = [NSString stringWithFormat:@"%@/IMG_BOLD_respcorr", work];
			//	[rindex saveAsKOImage:path];
				
				rindex = [rindex avgForLoop:bLp];
				p = [rindex data];
				for (i = 0; i < [avgLp dataLength]; i++) {
					if (p[i] < 0) p[i] = 0;
				}
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_rindex", work];
				[rindex saveAsKOImage:path];
				
				// ### calc weighted average (BOLD_avg2, BOLD_resp2)
				aimg = [bimg copy];
				[aimg multBy1dImage:rindex];
				aimg = [aimg avgForLoop:avgLp];
				[aimg subImage:[aimg avgForLoop:bLp]];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_avg2", work];
				[aimg saveAsKOImage:path];
				tmp_img = [aimg copy];
		//		[tmp_img gauss1DLP:0.4 forLoop:bLp];	// 0.08 BOLD
//		time_pca(tmp_img, 2);
//		exit(0);
				[tmp_img multBy1dImage:resp];
				tmp_img = [tmp_img avgForLoop:bLp];
				path = [NSString stringWithFormat:@"%@/IMG_BOLD_resp2", work];
				[tmp_img saveAsKOImage:path];

				// ### NCI filtering (3-5-2020)
				tmp_img = [aimg copy];
				[tmp_img gauss2DHP:48.0 / 64];
			path = [NSString stringWithFormat:@"%@/IMG_NCI_1", work];	// current ?
			[tmp_img saveAsKOImage:path];
			//	[tmp_img gauss1DLP:25.0 / 200 forLoop:bLp];
				[tmp_img gauss1DLP:25.0 / 200 forLoop:[tmp_img zLoop]];
				path = [NSString stringWithFormat:@"%@/IMG_NCI_avg", work];	// current ?
				[tmp_img saveAsKOImage:path];

		time_pca(tmp_img, 2);
		exit(0);

				[tmp_img accumFilt:10 forLoop:bLp];		// accum
			//	[tmp_img leakyIntWithTau:8.0 cutOff:0 forLoop:bLp];
				path = [NSString stringWithFormat:@"%@/IMG_NCI_acc", work];	// current ?
				[tmp_img saveAsKOImage:path];
				
	


				break;

			case 3 :

			// block avg (10 frames/sec, for 20sec)
				if (study <= 5) {
					bSize = 200;	// 800
				} else {
					bSize = 240;
				}
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];
				[img subImage:[img avgForLoop:[img zLoop]]];
				[img accumFilt:10 forLoop:[img zLoop]];		// accum

				




				break;

			case 6 :	// on-off phasor
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:tDim/2 chDim:2];
				for (i = 0; i < tDim; i++) {
					k = i / 100;
					if (k % 2 == 0) {	// on
						pp = [tmp_img data] + (i % 100 + k/2 * 100) * iDim;	// dest
					} else {			// off
						pp = [tmp_img data] + (i % 100 + (k-1)/2 * 100) * iDim + tDim/2 * iDim;	// dest
					}
					p  = [img data]		+ i * iDim;	// src
					// copy slice
					for (j = 0; j < iDim; j++) {
						pp[j] = p[j];
					}
				}
				break;
			}
			if (0) {
				path = [NSString stringWithFormat:@"%@/IMG_subsample", work];
				[img saveAsKOImage:path];

		//		time_ica(img, nComp);
				time_pca(img, 2);
				path = [NSString stringWithFormat:@"%@/IMG_subsample_pfilt", work];
				[img saveAsKOImage:path];
				var = [img varForLoop:[img zLoop]];
				path = [NSString stringWithFormat:@"%@/IMG_p_var", work];
				[var saveAsKOImage:path];
			}
    }   // @autoreleasepool

    return 0;
}

// add mask, path modifier
void
time_pca(RecImage *img, int nflt)
{
	int				tDim, iDim;
	int				xDim, yDim;
	int				i, j, k;
	RecImage		*img_a, *img_u;
	RecImage		*pca, *slice;
	float			*p, *q, *qq;
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

	// pca filter
	for (k = 0; k < nflt; k++) {
		// img -= pca(k) x s(k)
		slice = [pca sliceAtIndex:k forLoop:[pca zLoop]];
		p = [slice data];	// PCA
		q = [img_u data] + k * tDim;
		qq = [img data];
		for (i = 0; i < tDim; i++) {
			for (j = 0; j < iDim; j++) {
				qq[j] -= p[j] * q[i] * sres->s->data[k];
			}
			qq += xDim * yDim;
		}
	}
}

// ### started on 12-27-2018 (copied pca)
void
time_ica(RecImage *img, int nComp)
{
	int				tDim, iDim;
	int				xDim, yDim;
	int				i;
	RecImage		*img_a, *img_u;
	RecImage		*ica, *xw;
	NumMatrix		*A;
	NSDictionary	*res;

	tDim = [img zDim];
	xDim = [img xDim];
	yDim = [img yDim];
	iDim = xDim * yDim;
	img_a = [RecImage imageOfType:[img type] xDim:iDim yDim:tDim];
	ica = [RecImage imageOfType:[img type] xDim:xDim yDim:yDim zDim:nComp];
	
	[img_a copyImageData:img];
	[img_a dumpLoops];
	A = [img_a toMatrix];
//path = [NSString stringWithFormat:@"%@/IMG_A", work];
//[img_a saveAsKOImage:path];
//saveAsKOImage(A, path);

//	sres = Num_svd(A);
	res = [A icaForNC:nComp];

	img_a = [[res objectForKey:@"Vt"] toRecImage];
	[ica copyImageData:img_a];
	path = [NSString stringWithFormat:@"%@/IMG_PCA", work];
	[ica saveAsKOImage:path];

//	xw = [[[res objectForKey:@"U"] toRecImage] trans]; //Num_m_to_im(sres->U);
	xw = [[res objectForKey:@"U"] toRecImage];
	[xw trans];
	path = [NSString stringWithFormat:@"%@/IMG_U", work];
	[xw saveAsKOImage:path];
	[xw fft1d:[xw xLoop] direction:REC_INVERSE];
	path = [NSString stringWithFormat:@"%@/IMG_Uf", work];
	[xw saveAsKOImage:path];

	img_a = [[res objectForKey:@"Y"] toRecImage];
	[ica copyImageData:img_a];
	path = [NSString stringWithFormat:@"%@/IMG_ICA", work];
	[ica saveAsKOImage:path];

	xw = [(NumMatrix *)[res objectForKey:@"XW"] trans]; //Num_m_to_im(sres->U);
	path = [NSString stringWithFormat:@"%@/IMG_XW", work];
	[xw saveAsKOImage:path];


exit(0);
//	for (i = 0; i < 20; i++) {
//		printf("%d %f\n", i, sres->s->data[i]);
//	}
	tmp_img = [img_u copy];
	[tmp_img fft1d:[img_u xLoop] direction:REC_INVERSE];
	path = [NSString stringWithFormat:@"%@/IMG_U_ft", work];
	[tmp_img saveAsKOImage:path];
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

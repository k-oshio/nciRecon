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
//	1. LPF


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

void		time_pca(RecImage *img);
void		time_ica(RecImage *img, int nComp);
void		mk_kern(void);			// Gabor kernel
void		test_model(void);			// q & d model
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
    RecImage        *img, *tmp_img;
	RecImage		*mgmask, *spmask, *avg, *var;
	int				study, series;
	BOOL			bold;
	RecLoop			*tmLp, *phsLp, *avgLp, *outerLp, *blkLp, *xLp, *yLp;
	int				dda = 100;	// 10sec * 10phs
	int				xDim = 64, yDim = 64;
	int				tDim, iDim, bDim, oDim, nPixel;

	// PCA
	RecImage		*img_a, *img_u;
	Num_mat			*A;
	Num_svd_result	*sres;
	// PCA filter
	RecImage		*pca, *slice, *err;
	float			*p, *q, *pp, *qq;
	int				i, j, k, slc, blk;
	
    printf("nciRec12 (7-7-2018)\n");

    @autoreleasepool {
		img = nil;
	//	system("rm IMG_*");

// === test model
	if (0) {
		test_model();
		exit(0);
	}
// === read raw ===
		study = 6;
			// 1 : 2018-6-11-1
			// 2 : 2018-6-11-2
			// 3 : 2018-6-27
			// 4 : 2018-6-29
			// 5 : 2019-1-21
			// 6 : 2019-2-4
			// 7 : 2019-2-15
		series = 2;
			//==== 1.0 sec protocol
			// 1 : 10 * 10 dda, 10 * 380 (NCI)
			// 2 : 10 * 10 dda, 10 * 380 (BOLD with short stim) 
			// 3 : 10 * 10 dda, 10 * 380 (BOLD with continuous stim)
			//==== 1.2 sec protocol
			// 2,3: 10 * 10 dda, 12 * 320 
		bold = YES;
		avgMode = 5;	// 0: low-res, no block, 1: lo-res, 20 sec block, 2: hi-res, 20 sec block
						// 3: on/off phasor
						//						 5: lo-res, 24 sec block, 6: hi-res, 24 sec block
		nComp = 20;		// 

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
			break;	// study 6
		}
		work = [NSString stringWithFormat:@"%@/results/%d", base, series];
		path = [NSString stringWithFormat:@"%@/IMG_comb", work];
		img = [RecImage imageWithKOImage:path];

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
			if (i < iDim*36/64) {
				p[i] = 0;
			} else {
				p[i] = 1;
			}
		}
	// mag mask
		mgmask = [avg copy];
		[mgmask magnitude];
		[mgmask thresAt:0.2];
		[mgmask multByImage:spmask];
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
//		if (series == 2 || series == 3) {
		if (1) {
			// select sub-set for processing
			xLp = [RecLoop loopWithDataLength:xDim];
			yLp = [RecLoop loopWithDataLength:yDim];
			tDim = 3600; //1200;	// 3800 total, (100-on / 100-off) x 18
			iDim = xDim * yDim;
		//	avgLp = [img crop:avgLp to:380 startAt:dda];
			tmLp = [img crop:[img zLoop] to:tDim startAt:dda + 0];

			[img multByImage:mgmask];

			[img magnitude];
		//	[img phase];		// need to remove baseline phase first ...
			avg = [img avgForLoop:[img zLoop]];
			[img subImage:avg];
			path = [NSString stringWithFormat:@"%@/IMG_mg", work];
			[img saveAsKOImage:path];

		//	tmp_img = [img copy];
		//	[tmp_img gauss1DLP:0.7/5.0 forLoop:[img zLoop]];	// BOLD

		//	img = [tmp_img copy];	// BOLD
		//	[img subImage:tmp_img];	// NCI

		//	[img gauss1DHP:1.5/5.0 forLoop:[img zLoop] frac:1.0];	// NCI
		//	[img gauss1DBP:0.1 center:1.0/10.0 forLoop:[img zLoop]];	// NCI

			switch (avgMode) {
			case 0 :
			// time average (1 frame/sec)
				bSize = 10;
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img zLoop]];
				break;
			case 1 :
			// block avg (10 frames/sec, for 20sec)
				bSize = 200;	// 800
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];
				break;
			case 2 :
				bSize = 200;
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];

				[img gauss1DHP:0.5/5.0 forLoop:[img zLoop] frac:1.0];
				break;

			case 3 :	// on-off phasor
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
//	path = [NSString stringWithFormat:@"%@/IMG_subsample", work];
//	[tmp_img saveAsKOImage:path];
//	exit(0);
				break;

			case 5 :
			// block avg (10 frames/sec, for 24sec)
				printf("avg mode = %d\n", avgMode);
				bSize = 240;	// 800
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];
				break;
			case 6 :
				bSize = 240;
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];

				[img gauss1DHP:0.5/5.0 forLoop:[img zLoop] frac:1.0];
				break;


			}

//img = [img scaleZBy:4.0];
			path = [NSString stringWithFormat:@"%@/IMG_subsample", work];
			[img saveAsKOImage:path];
			tmp_img = [img oversample];
			path = [NSString stringWithFormat:@"%@/IMG_subsample_hi", work];
			[tmp_img saveAsKOImage:path];

	//		time_pca(img);
			time_ica(img, nComp);
	exit(0);

		// ============ phasor filter =================
			if (0) {
				mk_kern();
				img = phasor_filt(img, kc, ks, 11);
				path = [NSString stringWithFormat:@"%@/IMG_phasor", work];
				[img saveAsKOImage:path];

				avg = [img avgForLoop:[img zLoop]];
				path = [NSString stringWithFormat:@"%@/IMG_p_avg", work];
				[avg saveAsKOImage:path];
			//	var = [img varForLoop:[img zLoop] withMean:avg];
				var = [img varForLoop:[img zLoop]];
				path = [NSString stringWithFormat:@"%@/IMG_p_var", work];
				[var saveAsKOImage:path];

			// PCA
	//			time_pca(img);


				bDim = 100;
				tDim = [img zDim];
				phsLp = [RecLoop loopWithDataLength:nblk];
				avgLp = [RecLoop loopWithDataLength:tDim/nblk];
				tmp_img = [img copy];
				img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:phsLp, avgLp, yLp, xLp, nil];


				p = [tmp_img real];	// src
				q = [tmp_img imag];
				pp = [img real];	// dst
				qq = [img imag];

				for (i = 0; i < tDim; i++) {
					blk = i / bDim;
					k = i % bDim;
					for (j = 0; j < iDim; j++) {
						pp[(blk * bDim + k) * iDim + j] = p[i * iDim + j];
						qq[(blk * bDim + k) * iDim + j] = q[i * iDim + j];
					}
				}

			// phasor hist for max var pixel (32, 41)
			//	tmp_img = phasor_hist(img, 32, 41, 0.4);	// phasor peak (wd = 0.5)
				tmp_img = phasor_hist(img, 31, 44, 0.4);	// BOLD peak 1 (wd = 0.2)
			//	tmp_img = phasor_hist(img, 34, 43, 0.4);	// BOLD peak 2
			//	tmp_img = phasor_hist(img, 44, 35, 0.4);	// phasor peak 2
				[tmp_img gauss2DLP:0.15];

				path = [NSString stringWithFormat:@"%@/IMG_phasor_hist", work];
				[tmp_img saveAsKOImage:path];
			}
			
			tmp_img = [img copy];
			[tmp_img fft1d:[tmp_img zLoop] direction:REC_INVERSE];
			path = [NSString stringWithFormat:@"%@/IMG_fft", work];
			[tmp_img saveAsKOImage:path];

		// compare blocks visually (reorder)
			tmp_img = [RecImage imageOfType:RECIMAGE_REAL xDim:bDim yDim:tDim/bDim zDim:nPixel];
			qq = [tmp_img data];
			p = [img data];
			q = [mgmask data];
			slc = 0;
			for (k = 0; k < [mgmask dataLength]; k++) {
				if (q[k] == 0) continue;
				for (i = 0; i < tDim/bDim; i++) {
					for (j = 0; j < bDim; j++) {
						qq[slc * tDim + i * bDim + j] = p[(i * bDim + j) * 64*64 + k];
					}
				}
				slc++;
			}
			path = [NSString stringWithFormat:@"%@/IMG_blk", work];
			[tmp_img saveAsKOImage:path];
			
		exit(0);




		// block structure
			if (bold) {
				blkLp = [RecLoop loopWithDataLength:bDim];
				avgLp = [RecLoop loopWithDataLength:tDim/bDim];
				[img magnitude];
				[img multByImage:mgmask];
			//[img gauss1DLP:0.08 forLoop:[img zLoop]];	// =========
			// === block
				tmp_img = [img copy];
			//	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:avgLp, blkLp, yLp, xLp, nil];
			//	[img copyImageData:tmp_img];

			//	img = [img avgForLoop:avgLp];
				path = [NSString stringWithFormat:@"%@/IMG_blk", work];
				[img saveAsKOImage:path];
			} else {
			// baseline removal (necessary for unwrapping)
				tmp_img = [avg copy];
			//	[tmp_img gauss2DLP:0.1];
				[tmp_img conjugate];
				[tmp_img toUnitImage];
				[img multByImage:tmp_img];
				path = [NSString stringWithFormat:@"%@/IMG_baseline", work];
				[img saveAsKOImage:path];
			// take phase
				[img phase];
				[img multByImage:mgmask];
				if (0) {	// DC removal
					[img removeDC2dWithMask:mgmask];
				}
				if (0) {	// HPF
					[img gauss1DHP:0.05 forLoop:[img zLoop] frac:2.0];
				}
				if (1) {	// gamma filter (### not correct yet) -> make gamma filter method
					[img gammaFiltForLoop:[img zLoop] width:1.0];
				}
				path = [NSString stringWithFormat:@"%@/IMG_phase", work];
				[img saveAsKOImage:path];
			}

			tmp_img = [img copy];
			[tmp_img fft1d:[tmp_img zLoop] direction:REC_FORWARD];
			path = [NSString stringWithFormat:@"%@/IMG_fft", work];
			[tmp_img saveAsKOImage:path];
			
	//====== PCA
if (1) {
			tDim = [img zDim];
			img_a = [RecImage imageOfType:[img type] xDim:xDim * yDim yDim:tDim];
			pca = [RecImage imageOfType:[img type] xDim:xDim yDim:yDim zDim:tDim];
			
			[img_a copyImageData:img];
			A = Num_im_to_m(img_a);
			sres = Num_svd(A);
			img_a = Num_m_to_im(sres->Vt);
			[pca copyImageData:img_a];
			path = [NSString stringWithFormat:@"%@/IMG_PCA", work];
			[pca saveAsKOImage:path];
			img_u = Num_m_to_im(sres->U);
			[img_u trans];
			[img_u crop:[img_u yLoop] to:32 startAt:0];
			if (bold) {
			//	img_u = [img_u scale1dLoop:[img_u xLoop] to:[img_u xDim]/10];
			//	bDim /= 10;
			}
			path = [NSString stringWithFormat:@"%@/IMG_U", work];
			[img_u saveAsKOImage:path];
			for (k = 0; k < 20; k++) {
				printf("%d %f\n", k, sres->s->data[k]);
			}
			tmp_img = [img_u copy];
			[tmp_img fft1d:[img_u xLoop] direction:REC_INVERSE];
			path = [NSString stringWithFormat:@"%@/IMG_U_ft", work];
			[tmp_img saveAsKOImage:path];
}
		//====== PCA filter
			//====== BOLD
			if (bold) {
				// pca filter (select top 2)
				err = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim yDim:yDim zDim:bDim];
				for (k = 1; k < 2; k++) {
					// img -= temp_img(k) x s(k)
					slice = [pca sliceAtIndex:k forLoop:[tmp_img zLoop]];
					p = [slice data];	// PCA
					q = [img_u data] + k * bDim;
					qq = [err data];
					for (i = 0; i < bDim; i++) {
						for (j = 0; j < iDim; j++) {
							qq[j] += p[j] * q[i] * sres->s->data[k];
						}
						qq += iDim;
					}
				//	q  += tDim;
				}
				path = [NSString stringWithFormat:@"%@/IMG_PCAfilt", work];
				[err saveAsKOImage:path];
			} else {
			//====== NCI
				// img/err: input (IMG_phase), tDim/10:10:yDim:xDim
				// img_u:   U
				// pca: IMG_PCA, tDim:yDim:xDim
				err = [RecImage imageOfType:RECIMAGE_REAL withImage:img];
				tDim = [img_u xDim];
				for (k = 0; k < 1; k++) {
					// img -= temp_img(k) x s(k)
					slice = [pca sliceAtIndex:k forLoop:[tmp_img zLoop]];
					p = [slice data];
					q = [img_u data] + k * tDim;
					qq = [err data];
					for (i = 0; i < tDim; i++) {
						for (j = 0; j < iDim; j++) {
							qq[j] += p[j] * q[i] * sres->s->data[k];
						}
						qq += iDim;
					}
				//	q  += tDim;
					// chk self & var
				}
				[img subImage:err];
				path = [NSString stringWithFormat:@"%@/IMG_PCAfilt", work];
				[img saveAsKOImage:path];

			//==== block avg / var (block sub) ####
				phsLp = [RecLoop loopWithDataLength:10];
				avgLp = [RecLoop loopWithDataLength:10];
				blkLp = [RecLoop loopWithDataLength:2];
				tmp_img = [img copy];
				tDim = [img zDim];
				oDim = tDim / 10 / 10 / 2;
				printf("outerDim = %d\n", oDim);
				outerLp = [RecLoop loopWithDataLength:oDim];
				img = [RecImage imageOfType:RECIMAGE_REAL withLoops:outerLp, blkLp, avgLp, phsLp, yLp, xLp, nil];
				[img copyImageData:tmp_img];
				
				[img swapLoop:outerLp withLoop:blkLp];

				tmp_img = [img copy];
				avgLp = [RecLoop loopWithDataLength:tDim/2];
				img = [RecImage imageOfType:RECIMAGE_REAL withLoops:blkLp, avgLp, yLp, xLp, nil]; // wrong #### phsLp
				[img copyImageData:tmp_img];
				path = [NSString stringWithFormat:@"%@/IMG_blk", work];
				[img saveAsKOImage:path];

				avg = [img avgForLoop:avgLp];
				var = [img varForLoop:avgLp withMean:avg];

				path = [NSString stringWithFormat:@"%@/IMG_bavg", work];
				[avg saveAsKOImage:path];
				path = [NSString stringWithFormat:@"%@/IMG_bvar", work];
				[var saveAsKOImage:path];
				
		//		slice = [avg sliceAtIndex:0 forLoop:blkLp];
		//		tmp_img = [avg sliceAtIndex:1 forLoop:blkLp];
				slice = [var sliceAtIndex:0 forLoop:blkLp];
				tmp_img = [var sliceAtIndex:1 forLoop:blkLp];
				[slice subImage:tmp_img];
				[slice removeDC2dWithMask:mgmask];
				path = [NSString stringWithFormat:@"%@/IMG_bsub", work];
				[slice saveAsKOImage:path];

		//====== PCA2
				img = [slice copy];
				tDim = [img zDim];
				img_a = [RecImage imageOfType:[img type] xDim:xDim * yDim yDim:tDim];
				tmp_img = [RecImage imageOfType:[img type] xDim:xDim yDim:yDim zDim:tDim];
				
				[img_a copyImageData:img];
				A = Num_im_to_m(img_a);
				sres = Num_svd(A);
				img_a = Num_m_to_im(sres->Vt);
				[tmp_img copyImageData:img_a];
				path = [NSString stringWithFormat:@"%@/IMG_PCA2", work];
				[tmp_img saveAsKOImage:path];
				img_u = Num_m_to_im(sres->U);
				[img_u trans];
				[img_u crop:[img_u yLoop] to:32 startAt:0];
				if (bold) {
				//	img_u = [img_u scale1dLoop:[img_u xLoop] to:[img_u xDim]/10];
				//	bDim /= 10;
				}
				path = [NSString stringWithFormat:@"%@/IMG_U2", work];
				[img_u saveAsKOImage:path];
				for (k = 0; k < 10; k++) {
					printf("%d %f\n", k, sres->s->data[k]);
				}

			//==== k-filter
				tmp_img = [slice toDipole];
				path = [NSString stringWithFormat:@"%@/IMG_bsub_dp", work];
				[tmp_img saveAsKOImage:path];
			}
		}


// ===== straight
		if (series == 1) {
			xLp = [RecLoop loopWithDataLength:xDim];
			yLp = [RecLoop loopWithDataLength:yDim];
			tDim = 3600; //1200;	// 3800 total
			iDim = xDim * yDim;
			tmLp = [img crop:[img zLoop] to:tDim startAt:dda + tDim * 0];

			[img multByImage:mgmask];

			if (bold) {		// not "BOLD", but mag
				[img magnitude];
				tmp_img = [avg copy];
				[tmp_img magnitude];
				[img subImage:tmp_img];
				[img multByImage:mgmask];
			} else {		// phase
				// baseline removal (necessary for unwrapping)
				tmp_img = [avg copy];
				[tmp_img conjugate];
				[tmp_img toUnitImage];
				[img multByImage:tmp_img];
			//	path = [NSString stringWithFormat:@"%@/IMG_baseline", work];
			//	[img saveAsKOImage:path];
			// take phase
				[img phase];
				[img multByImage:mgmask];
				
			}
			// === time series filter ===
	//		[img gauss1DHP:0.2 forLoop:[img zLoop] frac:1.0];		// LPF (remove BOLD)
	//		[img gauss1DN:0.01 center:0.11 forLoop:[img zLoop]];	// notch (remove cardiac)
	//		[img gauss1DN:0.005 center:0.1 forLoop:[img zLoop]];	// notch (remove 1Hz)
//			[img t1DPF:100 forLoop:[img zLoop]];	// notch (remove cardiac)
			
			path = [NSString stringWithFormat:@"%@/IMG_in", work];
			[img saveAsKOImage:path];

			switch (avgMode) {
			case 0 :
			// time average (1 frame/sec)
				bSize = 10;
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:bSize chDim:[img zDim] / bSize];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img zLoop]];
				break;
			case 1 :
			// block avg (10 frames/sec, for 20sec)
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:200 chDim:[img zDim] / 200];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];
				break;
			case 2 :
				tmp_img = [RecImage imageOfType:[img type] xDim:[img xDim] yDim:[img yDim] zDim:10 chDim:[img zDim] / 10];
				[tmp_img copyImageData:img];
				img = [tmp_img avgForLoop:[tmp_img topLoop]];
				break;
			}

			path = [NSString stringWithFormat:@"%@/IMG_subsample", work];
			[img saveAsKOImage:path];
			

			time_pca(img);
			
		} // straight (series = 1)


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
	if (bold) {
	//	img_u = [img_u scale1dLoop:[img_u xLoop] to:[img_u xDim]/10];
	//	bDim /= 10;
	}
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


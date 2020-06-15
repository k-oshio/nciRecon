//
//	nci10 recon
//	12-15-2016		forked off from nciRec8
//	this is for Ogawa Lab experiments

//	=== plans ===
//	use dbg flag to switch dev / deploy
//	temporal freq filter
//	stim 0: no stim (after acq), 1: right before 90, 2: 100msec before 90


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

int
main(int ac, char *av[])
{
	LP_DIM			lp_dim;
    RecImage        *img, *img_t, *img_tf, *img_im;
	RecImage		*stim, *ref, *stim1, *stim2;
	RecImage		*pimg1, *pimg2;
	RecImage		*avg;
	RecImage		*mask, *hi;
	RecImage		*tmp;
//	RecLoop			*repLp;	// original loop structure in raw data
//	int				nRep;
    RecLoop         *avgLp, *phsLp, *slcLp;	// reorder to this
	RecLoop			*timeLp;
	RecLoop			*xLp, *yLp;
	int				nTime;
	int				nPhs;
	int				protocol;
	float			thres = 0.15;
	int				dda = 6;
	int				dbg = 1;	// dbg on : for tuning, cache on, read test data
								// dbg off: for actual processing, cache off, take argument as input file
	BOOL			useCache = NO;

    printf("nciRec10 (12-19-2016)\n");

    @autoreleasepool {
        system("rm IMG_*");
		img = nil;

// === read raw ===
		protocol = 0;	// decide using header ?
		img = [RecImage imageFromFile:@"IMGsav" relativePath:YES];
		if (img == nil || !useCache) {
			if (dbg) {
		// ### auditory ###
				img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0818a" andMeasData:@"meas_files/meas.out_0818a" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0818b" andMeasData:@"meas_files/meas.out_0818b" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0420a" andMeasData:@"meas_files/meas.out_0420a" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0420b" andMeasData:@"meas_files/meas.out_0420b" lpDim:&lp_dim];
		// ### visual ###
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_070627a" andMeasData:@"meas_files/meas.out_070627a" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_070627b" andMeasData:@"meas_files/meas.out_070627b" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_070919"  andMeasData:@"meas_files/meas.out_070919" lpDim:&lp_dim];
			//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_080109"  andMeasData:@"meas_files/meas.out_080109" lpDim:&lp_dim];
			} else {
			//	img = [RecImage imageWithMeasVD:@"meas.dat" lpDim:&lp_dim];
				img = [RecImage imageWithMeasAsc:@"meas.asc" andMeasData:@"meas.out" lpDim:&lp_dim];
			}
			// kx : nSamples    
			// ky : nLin
			// kz : nSlc
			// avg: nAcq
			// phs: nRep

			phsLp = [RecLoop findLoop:@"phs"];
			avgLp = [RecLoop findLoop:@"avg"];
			slcLp = [RecLoop findLoop:@"kz"];

			img = [img sliceAtIndex:0 forLoop:slcLp];

			[img swapLoop:phsLp withLoop:avgLp];
			[img dumpLoops];
[img saveAsKOImage:@"IMG_raw"];
			[img fft2d:REC_FORWARD];
			[img freqCrop];
			[img dumpLoops];
			[img saveAsKOImage:@"IMG_in"];
	// === corrections
			[img epiPcorr2];	// polynomial
		//	[img epiPcorr];	// 
			[img saveAsKOImage:@"IMG_epipcorr"];

//			img = [img rotShiftForEachOfLoop:slcLp scale:1.0];
		//	img = [img removeDdaAndMakePo2:dda];	// this changes result of rotShift
		//	[img rotShift];
		//	[img saveAsKOImage:@"IMG_rotShift"];
			[img baseline];
			[img saveAsKOImage:@"IMG_base"];
			[img saveToFile:@"IMGsav" relativePath:YES];	// cache this to save time
		}

// small FOV
		if (0) {
			RecLoop	*xLp, *yLp;
			int		siz = 32;	// 32
			int		xpos = 48;	// 48
			int		ypos = 72;	// 68
			xLp = [RecLoop loopWithDataLength:siz];
			yLp = [RecLoop loopWithDataLength:siz];
			[img replaceLoop:[img xLoop] withLoop:xLp offset:xpos];	// 32 : 60
			[img replaceLoop:[img yLoop] withLoop:yLp offset:ypos];
[img saveAsKOImage:@"IMG_sub.img"];
		}


// reorder into time series
//		img = [img removeDdaAndMakePo2:dda];	// no need for this... gaussHPF works for non-po2
		img_t = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:[img xDim] yDim:[img yDim]
					zDim:[avgLp dataLength] * [phsLp dataLength]];
		[img_t copyImageData:img];
		img = [img_t copy];
		
		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		
		[img saveAsKOImage:@"IMG_sav"];
		avg = [img avgForLoop:timeLp];
		[avg saveAsKOImage:@"IMG_avg"];
		hi = [avg oversample];
		[hi saveAsKOImage:@"IMG_avg_hi"];

// === mask
		mask = [img sdMaskForLoop:timeLp thres:thres];	// rect mask
		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
		[mask saveAsKOImage:@"IMG_mask"];

// === var (before time filter)
		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		img_im = [img copy];
		[img_im takeImagPart];
	//	[img_im takeRealPart];
	//	[img_im phase];
		[img_im saveAsKOImage:@"IMG_time_0"];

		tmp = [img_im sdForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_0"];
		if (0) {	// freq filter
			[img_im gauss1DHP:0.15 forLoop:timeLp frac:1.0];	// f0 = 0 (w = 0.25)
			[img_im saveAsKOImage:@"IMG_timefilt"];
		} else {	// PCA filter
			[img_im pcaFilt:20];	// currently this is just PCA. number of components is ignored
		}
exit(0);

// === var (after time filter)
		tmp = [img_im varForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_1"];
		[img_im saveAsKOImage:@"IMG_time_1"];

// === (2) kFilt
		img_t = [img_im kFilt:0];
		tmp = [img_t varForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_2"];
// === (3) apply mask
		[img_t saveAsKOImage:@"IMG_kft"];

		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
		[mask saveAsKOImage:@"IMG_mask"];

		[img_t multByImage:mask];
		[img_t saveAsKOImage:@"IMG_time"];

		tmp = [img_t varForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_mask"];

		img_tf = [img_t copy];
		[img_tf fft1d:timeLp direction:REC_INVERSE];
	//	[img_tf shift1d:timeLp];	// DC is 0
		[img_tf saveAsKOImage:@"IMG_tf"];

// == calc complex avg diff
		nPhs = 3;	// fixed ...
		// img_t cotains po2 images (without dda), phsArray is for full data set (including dda)
		ref   = [img_t phsSliceAtIndex:0 nPhs:nPhs dda:0];	// 256
		stim1 = [img_t phsSliceAtIndex:1 nPhs:nPhs dda:0];	// 128
		stim2 = [img_t phsSliceAtIndex:2 nPhs:nPhs dda:0];	// 128

// == calc p-image

		ref = [ref avgForLoop:[ref zLoop]];
		stim1 = [stim1 avgForLoop:[stim1 zLoop]];
		stim2 = [stim2 avgForLoop:[stim2 zLoop]];

		[ref   saveAsKOImage:@"IMG_ref"];
		[stim1 saveAsKOImage:@"IMG_stim1"];
		[stim2 saveAsKOImage:@"IMG_stim2"];

		stim1 = [stim1 subImage:ref];
		[stim1 saveAsKOImage:@"IMG_dif1"];
		stim2 = [stim2 subImage:ref];
		[stim2 saveAsKOImage:@"IMG_dif2"];

// == make combined result
		RecLoop	*stimLp = [RecLoop loopWithDataLength:2];
		xLp = [stim1 xLoop];
		yLp = [stim1 yLoop];
		stim = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:stimLp, yLp, xLp, nil];
		[stim copySlice:stim1 atIndex:0];
		[stim copySlice:stim2 atIndex:1];
		[stim saveAsKOImage:@"IMG_dif3"];

    }   // @autoreleasepool

    return 0;
}


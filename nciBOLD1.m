//
//	nciBOLD
//	12-23-2016		forked off from nciRec10

//	=== plans ===
//	for NIRS experiments (BOLD part)
//	dda already removed
//	10 off, 10 on, ... 0-1-0-1-0-1-0-1-0 (9 phases)
//	make loop structure, 4sec delay (2 samples, TR=2s), into 2 groups (40 on, 40 off)
//	*shift (linear), baseline (linear)



#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

int
main(int ac, char *av[])
{
	LP_DIM			lp_dim;
    RecImage        *img, *img_bold, *img_a;
	RecImage		*avg, *var, *hi, *pImg;
	RecImage		*tmp, *col;
    RecLoop         *avgLp, *phsLp, *slcLp;	// reorder to this
	RecLoop			*timeLp;
	int				nTime;
	int				xDim, yDim, zDim;

	int				protocol;
	int				evLen = 10;	// nOn == nOff
	int				nEvent = 4;
	int				nAvg;
	int				delay = 2;	// 4 sec
	int				dda = 0;

	int				dbg = 1;	// dbg on : for tuning, cache on, read test data
								// dbg off: for actual processing, cache off, take argument as input file

    printf("nciRec10 (12-22-2016)\n");

    @autoreleasepool {
        system("rm IMG_*");
		img = nil;

// === read raw ===
		protocol = 0;	// decide using header ?

		if (dbg) {
			img = [RecImage imageWithMeasVD:@"meas_files/meas_MID25_Checkerboard_block_FID37691.dat" lpDim:&lp_dim]; // 11/28
		//	img = [RecImage imageWithMeasVD:@"meas_files/meas_MID83_Checkerboard_block_FID35648.dat" lpDim:&lp_dim]; // 10/19
		} else {
			img = [RecImage imageWithMeasVD:@"meas.dat" lpDim:&lp_dim];
		}
		// kx : nSamples    
		// ky : nLin
		// kz : nSlc
		// avg: nAcq
		// phs: nRep

		phsLp = [RecLoop findLoop:@"phs"];
		avgLp = [RecLoop findLoop:@"avg"];
		slcLp = [RecLoop findLoop:@"kz"];

		[img reorderSlice]; // (Siemens)
		[img swapLoop:phsLp withLoop:slcLp];
		[img removePointLoops];
		[img dumpLoops];

		[img fft2d:REC_FORWARD];
		[img freqCrop];
		[img dumpLoops];
		[img saveAsKOImage:@"IMG_in"];
// === corrections
//		[img epiPcorr2];
		[img epiPcorrGE];
		[img saveAsKOImage:@"IMG_epipcorr"];

		[img magnitude];
		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		xDim = [img xDim];
		yDim = [img yDim];
		zDim = [[img topLoop] dataLength];

//		[img rotShift];	// doesn't work well -> drift corr
//		[img driftCorr];
		[img saveAsKOImage:@"IMG_base"];

		avg = [img avgForLoop:timeLp];
		[avg saveAsKOImage:@"IMG_avg"];
		hi = [avg oversample];
		[hi saveAsKOImage:@"IMG_avg_hi"];

// for ICA test
if (1) {
	Num_mat			*A;
	Num_svd_result	*sres;
	img = [img sliceAtIndex:9 forLoop:slcLp];
	[img saveAsKOImage:@"IMG_bold.sav"];
	img_a = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim * yDim yDim:nTime];
	[img_a copyImageData:img];
	A = Num_im_to_m(img_a);
	sres = Num_svd(A);
	saveAsKOImage(sres->U, @"IMG_U");
	img_a = Num_m_to_im(sres->Vt);
	[img copyImageData:img_a];
	[img saveAsKOImage:@"IMG_PCA"];
	Num_free_svd_result(sres);
	Num_free_mat(A);

	exit(0);
}
		
		nAvg = nEvent * evLen;
		avgLp = [RecLoop loopWithDataLength:nAvg];
		phsLp = [RecLoop loopWithDataLength:2];		// 2: on/off 
		// copy image data (reorder)
		img_bold = [img reorderBoldWithDda:dda nEvent:nEvent evLen:evLen slcLp:slcLp phsLp:phsLp avgLp:avgLp delay:delay];
		[img_bold saveAsKOImage:@"IMG_bold"];

		// difference of avg
		avg = [img_bold avgForLoop:avgLp];
		var = [img_bold varForLoop:avgLp withMean:avg];	// ### this is wrong... take var separately
		tmp = [avg difForLoop:phsLp stim:1 ref:0];
		[tmp saveAsKOImage:@"IMG_bold_dif"];

		// t-statistics (dif of avg)
		pImg = [img_bold pImageWithStim:1 ref:0 complex:NO];	// t or F-stat (for real/cpx)
		[pImg pvalThres:0.01];
		[pImg saveAsKOImage:@"IMG_pImg"];
		
		// fusion / color
		avg = [avg avgForLoop:phsLp];
		col = [avg fusionWithPimg:pImg gain:3000.0];
		col = [col oversample];
		[col saveAsKOImage:@"IMG_bold_fus"];

    }   // @autoreleasepool

    return 0;
}


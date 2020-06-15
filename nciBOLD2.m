//
//	nciBOLD2
//	3-12-2018		forked off from nciRecBOLD1
//					phase of GE EPI seq (BOLD set, or single stim set)
//	=== plans ===
//	just look at phase part of BOLD after unwrapping along time direction
//


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

int
main(int ac, char *av[])
{
	LP_DIM			lp_dim;
    RecImage        *img, *img_bold, *img_a, *img_d;
	RecImage		*avg, *var, *hi, *pImg, *mask;
	RecImage		*tmp, *col;
    RecLoop         *avgLp, *phsLp, *slcLp;	// reorder to this
	RecLoop			*timeLp;
	int				nTime;
	int				xDim, yDim, zDim;
	// PCA
	Num_mat			*A;
	Num_svd_result	*sres;

	int				dbg = 1;	// dbg on : for tuning, cache on, read test data
								// dbg off: for actual processing, cache off, take argument as input file

    printf("nciRecBOLD2 (3-12-2018)\n");

    @autoreleasepool {
        system("rm IMG_*");
		img = nil;

// === read raw ===
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

		[img pcorr1dForLoop:phsLp];
		[img saveAsKOImage:@"MG_timepcorr"];

	//	[img phase];
		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		xDim = [img xDim];
		yDim = [img yDim];
		zDim = [[img topLoop] dataLength];

//		[img rotShift];	// doesn't work well -> drift corr
//		[img driftCorr];
		[img saveAsKOImage:@"IMG_base"];

//		hi = [avg oversample];
//		[hi saveAsKOImage:@"IMG_avg_hi"];

// PCA
		img = [img sliceAtIndex:9 forLoop:slcLp];

		mask = [img sdMaskForLoop:[img zLoop] thres:0.1];	// rect mask
		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
		[mask saveAsKOImage:@"IMG_mask"];

		[img phase];
		[img multByImage:mask];

		[img saveAsKOImage:@"IMG_bold_p.sav"];
		img_a = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim * yDim yDim:nTime];
		[img_a copyImageData:img];
		A = Num_im_to_m(img_a);
		sres = Num_svd(A);
		saveAsKOImage(sres->U, @"IMG_U");
		img_a = Num_m_to_im(sres->Vt);
		[img copyImageData:img_a];
		[img saveAsKOImage:@"IMG_PCA"];

		img_d = [img kFilt:0];
		[img_d saveAsKOImage:@"IMG_kFilt"];

		Num_free_svd_result(sres);
		Num_free_mat(A);
		

    }   // @autoreleasepool

    return 0;
}


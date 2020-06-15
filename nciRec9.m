//
//	nci9 recon
//	10-19-2016		forked off from nciRec8
//	this is for NIRS experiments

//	=== plans ===
//	-> move test routines to single program (libtest.m ?)
//	-> test stat lib routines (F & t)	2-28-2017


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

void
delaySpectrum(int *p, int n)
{
	RecImage	*img;
	int			i;
	float		*pp;

	n = Rec_down2Po2(n);
	img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n];
	pp = [img data];
	for (i = 0; i < n; i++) {
		pp[i] = p[i];
	}
	[img fft1d:[img xLoop] direction:REC_INVERSE];
	[img saveAsKOImage:@"IMG_delay_f"];
}

// move to RecImageNCI.m when done ###
int *
read_protocol(NSString *path, int *npts, int *dda)
{
	NSArray			*lines;
    NSString		*content;
    NSError			*err;
	int				i, n, dd;
	float			*delay, del;
	int				phs, *phsArray;

    content = [NSString stringWithContentsOfFile:path encoding:NSASCIIStringEncoding error:&err];
    lines = [content componentsSeparatedByString:@"\n"];
	n = (int)[lines count];
	printf("%d lines found\n", n);

	delay = (float *)malloc(sizeof(float) * n);
	for (i = 0; i < n; i++) {
		delay[i] = [[lines objectAtIndex:i] floatValue];
	}
	dd = 0;
	for (i = 0; i < n; i++) {
		if (delay[i] < 0) {
			dd++;
		}
	}
	*dda = dd;
	n -= dd;

	phsArray = (int *)malloc(sizeof(int) * n);
	for (i = 0; i < n; i++) {
		del = delay[i + dd];
		if (del == 0) {
			phs = 0;
		} else
		if (del < 50.0) {
			phs = 1;
		} else {
			phs = 2;
		}
		phsArray[i] = phs;
printf("%d %d\n", i, phs);
	}
	*npts = n;
	free(delay);

	// spectrum of delay array
	delaySpectrum(phsArray, n);

	return phsArray;
}

int
main(int ac, char *av[])
{
	LP_DIM			lp_dim;	// ## slice loc ?
    RecImage        *img, *img_t, *img_im, *tmp, *img_tf, *avg;
//	RecImage		*pimg1, *pimg2;
	RecImage		*stim, *ref, *stim1, *stim2;
	RecImage		*mg;
	RecImage		*mask, *hi;
	RecLoop			*repLp;	// original loop structure in raw data
	int				nRep;
    RecLoop         *phsLp;	// reorder to this
	RecLoop			*xLp, *yLp;
	RecLoop			*timeLp;
	int				nPhs, nTime;
	int				protocol;
	int				*phsArray;
	float			thres = 0.15;
	int				dda;
	int				dbg = 1;	// dbg on : for tuning, cache on, read test data
								// dbg off: for actual processing, cache off, take argument as input file
//	int				statDbg = 0;
	int				timeFilt;
	BOOL			useCache = NO;

    printf("nciRec9 (10-19-2016)\n");

    @autoreleasepool {
        system("rm IMG_*");
		img = nil;


// === read raw ===
		protocol = 2;	// -> header
			// 1 : 10 dda, 260 x random(off, delay1 or 2)
			// 2 : 10 dda, 130 x (off, delay1, off, delay2)

	// protocol file contains lines like: -50:dda, float delay value in ms
		if (dbg) {
		// 10/19
			img = [RecImage imageWithMeasVD:@"meas_files/meas_MID89_DrOshio_fMRI_FID35654.dat" lpDim:&lp_dim];
			phsArray = read_protocol(@"protocol_35654.dat", &nRep, &dda);				
		// 11/28
		//	img = [RecImage imageWithMeasVD:@"meas_files/meas_MID32_DrOshio_fMRI_FID37698.dat" lpDim:&lp_dim];
		//	phsArray = read_protocol(@"protocol_37698.dat", &nRep, &dda);				
		} else {
			img = [RecImage imageWithMeasVD:@"meas.dat" lpDim:&lp_dim];
			phsArray = read_protocol(@"protocol.dat", &nRep, &dda);
				// same name for each case, placed within current dir
				// this file is created by processing csv source using readCsv
				// add nAvg, nPhs etc later
		}
		if (useCache) {
			img = [RecImage imageFromFile:@"IMGsav" relativePath:YES];
		} else {
			// kx : nSamples    
			// ky : nLin
			// kz : nSlc
			// avg: nAcq
			// phs: nRep

			[img fft2d:REC_FORWARD];
			[img freqCrop];
			[img dumpLoops];
			[img saveAsKOImage:@"IMG_in"];

	// === corrections
		//	[img epiPcorr2];
			[img epiPcorr];
			[img saveAsKOImage:@"IMG_epipcorr"];

//			img = [img rotShiftForEachOfLoop:slcLp scale:1.0];
			[img rotShift];
			[img saveAsKOImage:@"IMG_rotShift"];
			[img baseline];
			[img saveAsKOImage:@"IMG_base"];
			[img saveToFile:@"IMGsav" relativePath:YES];	// cache this to save time
		}
		[img saveAsKOImage:@"IMG_sav"];
		repLp = [RecLoop findLoop:@"phs"];
//		avg = [img avgForLoop:repLp];
//		[avg saveAsKOImage:@"IMG_avg"];
		hi = [avg oversample];
		[hi saveAsKOImage:@"IMG_avg_hi"];

	
// small FOV
		if (1) {
			RecLoop	*xLp, *yLp;
			int		siz = 32;	// 32
			int		xpos = 64;	// 48
			int		ypos = 105;	// 68
			xLp = [RecLoop loopWithDataLength:siz];
			yLp = [RecLoop loopWithDataLength:siz];
			[img replaceLoop:[img xLoop] withLoop:xLp offset:(xpos - siz/2)];	// 32 : 60
			[img replaceLoop:[img yLoop] withLoop:yLp offset:(ypos - siz/2)];
[img saveAsKOImage:@"IMG_sub.img"];
		avg = [img avgForLoop:repLp];
		[avg saveAsKOImage:@"IMG_avg"];
		hi = [avg oversample];
		[hi saveAsKOImage:@"IMG_avg_hi"];
		}

//exit(0);
	
// === mask 1 (mag mask, use phase variance)
		mask = [img sdMaskForLoop:repLp thres:thres];	// rect mask
		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
		[mask saveAsKOImage:@"IMG_mask_1"];

// === make loop structure
//	input:	nAcq = 1,	nRep = 530	-> time series
//	output: nAvg = 128, nPhs = 4	-> inter-group difference

//		dda = 10;	// read from protocol file
		nPhs = 4;
		phsLp = [RecLoop loopWithDataLength:nPhs];

// === test tag (update with delay array)
//	img contains all data (including dda) at this point
//	dda is set by read_protocol()
		if (0) {
			[img addTestTag:dda nPhs:nPhs phsArray:phsArray amp:100];	// ### not done yet
			[img saveAsKOImage:@"img_tag.img"];
		}

// === time series / base image
		img = [img removeDdaAndMakePo2:dda];

// ===== time series filters ... order is significant =======

// === (1) HPF in time direction
		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		img_im = [img copy];
		[img_im takeImagPart];
		[img_im saveAsKOImage:@"IMG_time_0"];

		tmp = [img_im sdForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_0"];

		timeFilt = 4;	// 1 seems to be ok...
		switch (timeFilt) {
		case 0 :	// skip
		default :
			break;
		case 1 :
			[img_im gauss1DHP:0.15 forLoop:timeLp frac:1.0];	// f0 = 0 (w = 0.25)
			break;
		case 2 :
			[img_im cos1dForLoop:timeLp cyc:2 power:1 lowPass:NO];		// or f0 = 0 & nyq
			break;
		case 3 :
			[img_im cos1dForLoop:timeLp cyc:4 power:4 lowPass:YES];
			[img_im gauss1DHP:0.15 forLoop:timeLp frac:1.0];	// f0 = 0 (w = 0.25)
			break;
		case 4 :
			[img_im pcaFilt:256];	// PCA type variance filter
			break;
		}
		[img_im saveAsKOImage:@"IMG_time_1"];

		tmp = [img_im sdForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_1"];

exit(0);

// === (2) kFilt
		img_t = [img_im kFilt:0];
		[img_t saveAsKOImage:@"IMG_time_2"];

		tmp = [img_t sdForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_2"];

// === (3) apply mask

//		mask = [img sdMaskForLoop:repLp thres:thres];	// rect mask
//		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
		[mask saveAsKOImage:@"IMG_mask"];

		[img_t multByImage:mask];
		[img_t saveAsKOImage:@"IMG_time"];
		tmp = [img_t varForLoop:timeLp];
		[tmp saveAsKOImage:@"IMG_var_mask"];

// === freq
		img_tf = [img_t copy];
		[img_tf fft1d:timeLp direction:REC_INVERSE];
		[img_tf shift1d:timeLp];	// DC is 0
		[img_tf saveAsKOImage:@"IMG_tf"];

// === correlation in time
		if (0) {
			tmp = [img_tf copy];
			[tmp magnitudeSq];
			[tmp fft1d:timeLp direction:REC_FORWARD];
			mg = [tmp sliceAtIndex:nTime/2 forLoop:timeLp];
			[tmp divImage:mg];
			[tmp saveAsKOImage:@"IMG_autoc"];
		}

// == calc complex avg diff
		// img_t cotains po2 images (without dda), phsArray is for full data set (including dda)
		ref   = [img_t phsSliceAtIndex:0 withPhsArray:phsArray dda:dda];	// 256
		stim1 = [img_t phsSliceAtIndex:1 withPhsArray:phsArray dda:dda];	// 128
		stim2 = [img_t phsSliceAtIndex:2 withPhsArray:phsArray dda:dda];	// 128

// == diff of mean
		ref   = [ref   avgForLoop:[ref   zLoop]];
		stim1 = [stim1 avgForLoop:[stim1 zLoop]];
		stim2 = [stim2 avgForLoop:[stim2 zLoop]];
		[stim1 subImage:ref];
		[stim2 subImage:ref];
		[stim1 saveAsKOImage:@"IMG_dif1"];
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


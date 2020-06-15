//
//	nciRec11 recon
//	6-12-2018		forked off from nciRec9
//	this is for NIRS experiments
//	new scheme using GR-EPI and long delay
//	first phase only (second phase is done in nciRec12)

//	=== plans ===


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

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

	return phsArray;
}

int
main(int ac, char *av[])
{
	LP_DIM			lp_dim;	// ## slice loc ?
    RecImage        *img, *avg, *prof;
	RecLoop			*repLp;	// original loop structure in raw data
	RecLoop			*chLp;
	int				study, series;

	// locations
	NSString		*base, *path, *work;

    printf("nciRec11 (6-12-2018)\n");

    @autoreleasepool {
		img = nil;
		system("rm IMG_*");

// === read raw ===
		study = 2;
			// 1 : 2018-6-11-1
			// 2 : 2018-6-11-2
			// 3 : 2018-6-27 (rot)
			// 4 : 2018-6-29 (rot)
			// 5 : 2019-1-21 ser 1 bad data
			// 6 : 2019-2-4
			// 7 : 2019-2-15
		series = 1;
			// 1 : 10 * 10 dda, 10 * 380 (NCI)
			// 2 : 10 * 10 dda, 10 * 380 (BOLD with short stim) 
			// 3 : 10 * 10 dda, 10 * 380 (BOLD with continuous stim)
			// ===== 1.2 sec protocol
			// 2,3 : 10 * 10 dda, 12 * 320

	//	if (useCache) {
	//		img = [RecImage imageFromFile:@"IMGsav" relativePath:YES];
	//	} else {
		// protocol file contains lines like: -50:dda, float delay value in ms
printf("read raw\n");
		switch (study) {
		case 1:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-1";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00032_FID06871_Checkerboard_block_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00034_FID06873_Checkerboard_block_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00036_FID06875_Checkerboard_block_S3.dat", base];
				break;
			}
			break;	// study 1
		case 2:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0611-2";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00067_FID06906_Checkerboard_block_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00069_FID06908_Checkerboard_block_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00071_FID06910_Checkerboard_block_S3.dat", base];
				break;
			}
			break;	// study 2
		case 3:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0627";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00049_FID07509_Checkerboard_block_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00051_FID07511_Checkerboard_block_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00055_FID07515_Checkerboard_block_S3.dat", base];
				break;
			}
			break;	// study 3
		case 4:
			base = @"/Users/oshio/images/NCI/NIRS/2018-0629";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00030_FID07768_Checkerboard_block_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00032_FID07770_Checkerboard_block_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00034_FID07772_Checkerboard_block_S3.dat", base];
				break;
			}
			break;	// study 4
		case 5:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0121";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00044_FID22723_Checkerboard_block(1_2s)_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00066_FID22745_Checkerboard_block(1_2s)_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00068_FID22747_Checkerboard_block(1_2s)_S1.dat", base];
				break;
			}
			break;	// study 5
		case 6:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0204";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00125_FID23778_Checkerboard_block(1_2s)_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00127_FID23780_Checkerboard_block(1_2s)_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00129_FID23782_Checkerboard_block(1_2s)_S3_state_rip.dat", base];
				break;
			}
			break;	// study 6
		case 7:
			base = @"/Users/oshio/images/NCI/NIRS/2019-0215";
			switch (series) {
			case 1:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00050_FID24646_Checkerboard_block(1_2s)_S1.dat", base];
				break;
			case 2:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00052_FID24648_Checkerboard_block(1_2s)_S2.dat", base];
				break;
			case 3:
				path = [NSString stringWithFormat:@"%@/data/meas_MID00054_FID24650_Checkerboard_block(1_2s)_S3.dat", base];
				break;
			}
			break;	// study 7
		}
		work = [NSString stringWithFormat:@"%@/results/%d", base, series];

		if (0) {
			[RecImage readMetaVD:path lpDim:&lp_dim];
			exit(0);
		}

		img = [RecImage imageWithMeasVD:path lpDim:&lp_dim];
		[img removePointLoops];
		[img dumpLoops];
//		path = [NSString stringWithFormat:@"%@/IMG_raw", work];
//		[img saveAsKOImage:path];

//		path = [NSString stringWithFormat:@"%@/Chckerboard_event_250ms_sparseTrigger_Cum-06112-1_TimingLog.txt", base];
//		phsArray = read_protocol(path, &nRep, &dda);
			// same name for each case, placed within current dir
			// this file is created by processing csv source using readCsv
			// add nAvg, nPhs etc later

		// kx : nSamples    
		// ky : nLin
		// kz : nSlc
		// avg: nAcq
		// phs: nRep

printf("fft\n");
		[img fft2d:REC_FORWARD];
		[img freqCrop];
		[img dumpLoops];
	//	path = [NSString stringWithFormat:@"%@/IMG_in", work];
	//	[img saveAsKOImage:path];

// === corrections
printf("epiPcorr3\n");
		[img epiPcorr3];
//		path = [NSString stringWithFormat:@"%@/IMG_epipcorr", work];
//		[img saveAsKOImage:path];

//		[img saveAsKOImage:@"IMG_sav"];
		repLp = [RecLoop findLoop:@"rep"];
		chLp = [RecLoop findLoop:@"ch"];

printf("coil prof calc\n");
		prof = [img coilProfile2DForLoop:chLp];
printf("combine\n");
		img = [img combineForLoop:chLp withProfile:prof];
		if (study == 3 || study == 4) {
			[img rotate:1];
		}
		path = [NSString stringWithFormat:@"%@/IMG_comb", work];
		[img saveAsKOImage:path];

printf("avg\n");
		avg = [img avgForLoop:repLp];
		path = [NSString stringWithFormat:@"%@/IMG_avg", work];
		[avg saveAsKOImage:path];

    }   // @autoreleasepool

// === make separate prog for actual processing ===

    return 0;
}


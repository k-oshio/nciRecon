//
//	nci8 recon
//	phase-based version
//	4-13-2016	K. Oshio
//	--> now current
//	5-5-2016	first definite result !!! -> X
//	6-25-2016	phase variance mask -> result -> X
//	10-5-2016	restarted... half gaussian mask works
//	10-8-2016	nphs = 2 for freq analysis

//	=== plans ===
//	freq analysis -> move to nci10


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

int
main()
{
    RecImage        *img, *tmp, *avg, *pimg, *col;
    RecImage        *ref, *stim, *phs, *mask;
    RecLoop         *avgLp, *phsLp, *slcLp;
	RecLoop			*blkLp, *timeLp, *ftLp;
	RecLoopControl	*lc;
    LP_DIM          lp_dim;
	int				blk = 2;	// 3
	int				moment = 0;	// 0:dipole, 1:quadrupole
	BOOL			useImag = NO;
	BOOL			oldBaseline = NO;
	BOOL			useCache = NO;
	BOOL			rot;
	BOOL			rotFirst = NO;
	BOOL			Po2 = YES;
	float			thres;

    printf("nciRec8 (6-29-2016)\n");

	thres = 0.12;	// 0.15
	rot = YES;

    @autoreleasepool {
        system("rm IMG_*");
		img = nil;
		if (useCache) {
			img = [RecImage imageFromFile:@"IMGsav" relativePath:YES];
			[img dumpLoops];
		}

// === first step ===
		if (img == nil) {
		// read raw / ft / phase correction / registration
		//  img = [RecImage imageWithMeasAsc:@"meas.asc" andMeasData:@"meas.out" lpDim:&lp_dim];
		// ### auditory ###
		// #  img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0511" andMeasData:@"meas_files/meas.out_0511" lpDim:&lp_dim];
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0818a" andMeasData:@"meas_files/meas.out_0818a" lpDim:&lp_dim];
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0818b" andMeasData:@"meas_files/meas.out_0818b" lpDim:&lp_dim];
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0420a" andMeasData:@"meas_files/meas.out_0420a" lpDim:&lp_dim];
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0420b" andMeasData:@"meas_files/meas.out_0420b" lpDim:&lp_dim];
		// ### visual ###
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_nci7" andMeasData:@"meas_files/meas.out_nci7" lpDim:&lp_dim];
		//	img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_071010a" andMeasData:@"meas_files/meas.out_071010a" lpDim:&lp_dim];
			img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_070919" andMeasData:@"meas_files/meas.out_070919" lpDim:&lp_dim];
		    // kx : nSamples    
			// ky : nLin
			// kz : nSlc
			// avg: nAcq
			// phs: nRep
	//printf("lp_dim.repIsOuter = %d\n", lp_dim.repIsOuter);
			[img fft2d:REC_FORWARD];
			[img freqCrop];
			[img dumpLoops];
			[img saveAsKOImage:@"IMG_in"];

			avgLp = [RecLoop findLoop:@"avg"];
			phsLp = [RecLoop findLoop:@"phs"];
			slcLp = [RecLoop findLoop:@"kz"];

	// crop to Power of 2 (phs & avg)
			if (Po2) {
				avgLp = [img cropToPo2:avgLp];
				[RecLoop replaceLoopNamed:@"avg" withLoop:avgLp];
				phsLp = [img crop:phsLp to:2];
				[RecLoop replaceLoopNamed:@"phs" withLoop:phsLp];
			//	[img saveAsKOImage:@"IMG_po2"];
				[img dumpLoops];
			//	exit(0);
			}

// dbg ... set marker
{
	float	*p = [img data] + [img dataLength];
	int		i, n = [avgLp dataLength];
	int		skip = [img skipSizeForLoop:avgLp];

	p += 8256;	// (4, 4)
	for (i = 0; i < n; i++) {
		*p += 1000;
		p += skip;
	}
}
//[img saveAsKOImage:@"IMG_in"];
//exit(0);

		[img swapLoop:phsLp withLoop:avgLp];
//		[img dumpLoops];

			// old (ref)
			if (oldBaseline) {	// for comparison with old method
				[img epiPcorr];
				[img saveAsKOImage:@"IMG_epiPcorr"];
				[img baseline];
			} else {
			// new baseline phase removal
				[img epiPcorr2];
				[img saveAsKOImage:@"IMG_epipcorr"];
				if (rotFirst) {
					if (rot) {
						printf("rotShift\n");
						img = [img rotShiftForEachOfLoop:slcLp scale:1.0];
						[img saveAsKOImage:@"IMG_rotShift"];
					}
				}
				[img baseline2Avg:avgLp phs:phsLp];
			}
			[img rotImage:lp_dim.inplaneRot];   // Siemens
			[img reorderSlice];                 // Siemens

	        [img saveAsKOImage:@"IMG_base"];

//			[img dumpPhaseSDForLoop:avgLp image:@"IMG_PHSD" histo:@"hist.txt"];

			if (!rotFirst) {
				if (rot) {
					printf("rotShift\n");
					img = [img rotShiftForEachOfLoop:slcLp scale:1.0];
					[img saveAsKOImage:@"IMG_rotShift"];
				}
			}
			[img saveToFile:@"IMGsav" relativePath:YES];	// cache this to save time
		}

		avg = [img avgForLoop:avgLp];
		[avg scaleToVal:1000.0];		// scale for fusion
		[avg saveAsKOImage:@"IMG_amp_mn"];
		tmp = [img sdForLoop:avgLp withMean:avg];
		[tmp saveAsKOImage:@"IMG_amp_var"];

		mask = [img sdMaskForLoop:avgLp andLoop:phsLp thres:thres];	// rect mask
		[mask hGauss2DLP:0.5];	// "Half" gaussian filter

// ### dbg
	if ([mask minVal] < 0) {
		printf(" ##### neg value detected\n");
	}
		[mask saveAsKOImage:@"IMG_mask"];

		phs = [img copy];
		if (useImag) {	// imag part
			[phs takeImagPart];
			[phs scaleByImage:mask];
		} else {	// phase
			[phs phase];
			[phs scaleByImage:mask];	// float mask
		}
// [phs saveAsKOImage:@"IMG_phs_m1"];

	// rho-direction filter
		[phs gauss2DHP:0.2 frac:1.0];
//		[phs gauss2DLP:0.6];		// 0.6
// [phs saveAsKOImage:@"IMG_phs_m2"];
		tmp = [phs avgForLoop:avgLp];
		[tmp saveAsKOImage:@"IMG_phs_mn"];
		tmp = [phs sdForLoop:avgLp withMean:tmp];
		[tmp saveAsKOImage:@"IMG_phs_var"];

	// spatial FT
		tmp = [phs copy];
		[tmp fft2d:REC_INVERSE];
		tmp = [tmp avgForLoop:avgLp];
		[tmp saveAsKOImage:@"IMG_phs_k"];

	// time FT
		lc = [phs control];
		timeLp = [lc combineLoop:avgLp andLoop:phsLp];
		tmp = [RecImage imageOfType:[phs type] withControl:lc];
		[tmp copyLoopsOf:phs];
		[tmp copyImageData:phs];
		[tmp saveAsKOImage:@"IMG_phs_time"];
//		ftLp = [tmp crop:timeLp to:128];
		ftLp = [tmp zeroFillToPo2:timeLp];
		[tmp fft1d:ftLp direction:REC_INVERSE];
		[tmp saveAsKOImage:@"IMG_phs_freq"];

	// time filter
		lc = [phs control];
		timeLp = [lc combineLoop:avgLp andLoop:phsLp];
		tmp = [RecImage imageOfType:[phs type] withControl:lc];
		[tmp copyImageData:phs];

		[tmp gauss1DHP:0.5 forLoop:timeLp frac:1.0];
		if (Po2) {
			[tmp gauss1DBP:0.04 center:0.5 forLoop:timeLp];
		} else {
			[tmp gauss1DBP:0.04 center:0.333 forLoop:timeLp];
		}

		[tmp saveAsKOImage:@"IMG_phs_flt"];
		ftLp = [tmp zeroFillToPo2:timeLp];
		[tmp fft1d:ftLp direction:REC_INVERSE];
		[tmp saveAsKOImage:@"IMG_phs_fltfreq"];
		tmp = [tmp sdForLoop:ftLp];
		[tmp saveAsKOImage:@"IMG_phs_fltfreq_sd"];
		
		[phs copyImageData:tmp];
		
exit(0);


	// without kFilt (real mean dif, real var dif)
		//======
        ref = [phs sliceAtIndex:0 forLoop:phsLp];       // first phs
     //   stim = [phs removeSliceAtIndex:0 forLoop:phsLp];
		stim = [phs copy];
		[stim removeSliceAtIndex:0 forLoop:phsLp];
	
       // pimg = [stim pImageWithRef:ref forLoop:avgLp thres:0.05 mode:NCI_REAL]; //ttest:YES];	// 0.1
		// [pimg maskWithImage:mask];

		pimg = [stim meanDiffWithRef:ref forLoop:avgLp];
		[pimg saveAsKOImage:@"IMG_mnDif0"];
		pimg = [stim sdDiffWithRef:ref forLoop:avgLp];
		[pimg saveAsKOImage:@"IMG_sdDif0"];
if (0) {
		//======
		col = [avg fusionWith:pimg gain:30.0 mode:0x01];
		[col saveAsKOImage:@"IMG_fus0"];
        col = [col oversample];
        [col saveAsKOImage:@"IMG_fus0_hi"];
}
// ### var dif


	// k-filter (th-direction)
        phs = [phs kFilt:moment];
		[phs saveAsKOImage:@"IMG_moment"];
		tmp = [phs varForLoop:avgLp];
		[tmp saveAsKOImage:@"IMG_phs_cvar"];

		//======
		tmp = [RecImage imageOfType:[phs type] withControl:lc];
		[tmp copyLoopsOf:phs];
		[tmp copyImageData:phs];
		[tmp saveAsKOImage:@"IMG_phs_ctime"];
		ftLp = [tmp crop:timeLp to:256];
		[tmp fft1d:ftLp direction:REC_INVERSE];
		[tmp saveAsKOImage:@"IMG_phs_cfreq"];

		//======
        ref = [phs sliceAtIndex:0 forLoop:phsLp];       // first phs
    //    stim = [phs removeSliceAtIndex:0 forLoop:phsLp];     // 2 and 3 (preserve phs for blk) ###
		stim = [phs copy];
		[stim removeSliceAtIndex:0 forLoop:phsLp];

		pimg = [stim meanDiffWithRef:ref forLoop:avgLp];
		[pimg saveAsKOImage:@"IMG_mnDif1"];
		pimg = [stim sdDiffWithRef:ref forLoop:avgLp];
		[pimg saveAsKOImage:@"IMG_sdDif1"];

        pimg = [stim pImageWithRef:ref forLoop:avgLp thres:0.05  mode:NCI_CPX_MEAN];	// 0.1
		[pimg maskWithImage:mask];
		//======
        [pimg saveAsKOImage:@"IMG_pimg_t"];

		col = [avg fusionWith:pimg gain:100.0 mode:0x01];
		[col saveAsKOImage:@"IMG_fus2"];
        col = [col oversample];
        [col saveAsKOImage:@"IMG_fus2_hi"];
exit(0);
    // ===== var-test =====
        pimg = [stim pImageWithRef:ref forLoop:avgLp thres:0.05  mode:NCI_CPX_VAR];
		[pimg maskWithImage:mask];
        [pimg saveAsKOImage:@"IMG_pimg_f"];
        col = [avg fusionWith:pimg gain:50.0 mode:0x01];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
        col = [col oversample];
        [col saveAsKOImage:@"IMG_fus4_hi"];
    // ===== blk avg (merge with mean-test when done) =====
        if (blk > 0 && [avgLp dataLength] > 20) {
			// create blk loop
			// breakup avg -> blk x avg
			blkLp = [phs breakLoop:avgLp blockSize:blk dummy:10];
			//======
			ref = [phs sliceAtIndex:0 forLoop:phsLp];       // first phs
		//	stim = [phs removeSliceAtIndex:0 forLoop:phsLp];     // 2 and 3 (preserve phs for blk)
			stim = [phs copy];
			[stim removeSliceAtIndex:0 forLoop:phsLp];
			pimg = [stim pImageWithRef:ref forLoop:blkLp thres:0.05  mode:NCI_CPX_MEAN];
			[pimg maskWithImage:mask];
			//======
			[pimg saveAsKOImage:@"IMG_blk_pimg"];
			col = [avg fusionWith:pimg gain:100.00 * sqrt(blk) mode:0x01];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
			[col saveAsKOImage:@"IMG_fus3"];
			col = [col oversample];
			[col saveAsKOImage:@"IMG_fus3_hi"];
        }
        avg = [avg oversample];
        [avg saveAsKOImage:@"IMG_avg_hi"];

    }   // @autoreleasepool

    return 0;
}


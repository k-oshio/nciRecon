//
//	nci7 recon
//	-> reproduce functionarity of C version nci_nci7_sc.c
//

#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

int
main()
{
    RecImage        *img, *avg, *var, *pimg, *col, *psum;
    RecImage        *ref, *phs;
    RecLoop         *avgLp, *phsLp, *slcLp, *blkLp;
    LP_DIM          lp_dim;
    float           wd = 2.0;
    int             blk = 0;    // 8
	int				moment = 1;	// 0:dipole, 1:quadrupole
	int				thres = 1; //2;	// 0: 0.05, 1: 0.01, 2: 0.005

    printf("nciRec7 (1-16-2016)\n");
	printf("moment:%d, thres:%d\n", moment, thres);

    @autoreleasepool {
        system("rm IMG_*");
    // read raw / ft / phase correction / registration
    //    img = [RecImage imageWithMeasAsc:@"meas.asc" andMeasData:@"meas.out" lpDim:&lp_dim];
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0511" andMeasData:@"meas_files/meas.out_0511" lpDim:&lp_dim];
        img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_nci2" andMeasData:@"meas_files/meas.out_nci2" lpDim:&lp_dim];
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0818b" andMeasData:@"meas_files/meas.out_0818b" lpDim:&lp_dim];
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_nci7" andMeasData:@"meas_files/meas.out_nci7" lpDim:&lp_dim];
       // kx : nSamples    
        // ky : nLin
        // kz : nSlc (point, removed)
        // avg: nAcq
        // phs: nRep
        avgLp = [RecLoop findLoop:@"avg"];
        phsLp = [RecLoop findLoop:@"phs"];
        slcLp = [RecLoop findLoop:@"kz"];
        [img fft2d:REC_FORWARD];
        [img freqCrop];
        [img epiPcorr];
//    [img pcorr];
        [img baseline];
        [img rotImage:lp_dim.inplaneRot];   // Siemens
        [img reorderSlice];                 // Siemens

        [img saveAsKOImage:@"IMG_base"];
        printf("rotShift\n");
        img = [img rotShiftForEachOfLoop:slcLp scale:1.0];
        [img saveAsKOImage:@"IMG_rotShift"];

        avg = [img avgForLoop:avgLp];       // avg for "avg" loop
        avg = [avg avgForLoop:phsLp];
        [avg saveAsKOImage:@"IMG_avg"];    // straight recon [phs avg kz ky kx]
        var = [img varForLoop:avgLp];
        var = [var avgForLoop:phsLp];
        [var sqroot];
        [var saveAsKOImage:@"IMG_var"];

    //    col = [avg fusionWith:var gain:3.0 mode:0x00];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
    //    col = [col oversample];
     //   [col saveAsKOImage:@"IMG_var"];
    // remove real part
        [img removeRealPart];
    
//        if (blk > 0) {    // blk avg (not corrent yet)
//            avgLp = [img breakLoop:avgLp blockSize:8 dummyShot:1];
//            [img saveAsKOImage:@"IMG_blk"];
//        }

    // k-filter
        ref = [img sliceAtIndex:0 forLoop:phsLp];       // first phs
        ref = [ref kFilt:moment];
        [ref gauss2DLP:wd]; // smoothing
[ref saveAsKOImage:@"IMG_ref"];

   //    phs = [img removeSliceAtIndex:0 forLoop:phsLp];     // 2 and 3
		phs = [img copy];
		[phs removeSliceAtIndex:0 forLoop:phsLp];	// phs has new (modified) phase loop
        phs = [phs kFilt:moment];
        [phs gauss2DLP:wd]; // smoothing
[phs saveAsKOImage:@"IMG_phs"];

        img = [img kFilt:moment];
        [img swapLoop:phsLp withLoop:slcLp];	// old (real) phsLp
        [img swapLoop:phsLp withLoop:avgLp];
        [img saveAsKOImage:@"IMG_time"];

    // ===== mean-test =====
//        pimg = [phs pImageWithRef:ref forLoop:avgLp thres:thres ttest:YES];
//        [pimg saveAsKOImage:@"IMG_pimg_t"];
//        col = [avg fusionWith:pimg gain:0.03 mode:0x01];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
//        col = [col oversample];
//        [col saveAsKOImage:@"IMG_fus2_hi"];

    // ===== var-test =====
//        pimg = [phs pImageWithRef:ref forLoop:avgLp thres:thres ttest:NO];   // mode:4
//        [pimg saveAsKOImage:@"IMG_pimg_f"];
//        col = [avg fusionWith:pimg gain:0.03 mode:0x01];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
//        col = [col oversample];
//        [col saveAsKOImage:@"IMG_fus4_hi"];

        if (blk > 0) {
        //    pimg = [pimg oversample];
        //    [pimg saveAsKOImage:@"IMG_pimg_hi"];
            blkLp = [RecLoop findLoop:@"blk"];
            psum = [pimg avgForLoop:blkLp];
        //    psum = [pimg tip:5 forLoop:blkLp];
        //    psum = [pimg mipForLoop:blkLp];
            [psum saveAsKOImage:@"IMG_psum"];
            col = [avg fusionWith:psum gain:0.03 mode:0x11];
            col = [col oversample];
            [col saveAsKOImage:@"IMG_fus5_hi"];
        }

        avg = [avg oversample];
        [avg saveAsKOImage:@"IMG_avg_hi"];

    }   // @autoreleasepool
    return 0;
}


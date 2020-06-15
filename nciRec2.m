//
//	nci2 recon
//  auditory, SE-EPI
//  last updated on 2-14-2014
//

// === protocol ===
// unified version. handles auditory / se-epi / ??
//  (chk error with 03-0508a)
//      

#import <RecKit/RecKit.h>
#import "RecImageNCI.h"

int
main(int ac, char *av[])
{
    NSString        *asc_path, *out_path;
    RecImage        *img, *avg, *pimg, *col;
    RecImage        *ref, *phs;
    RecLoop         *avgLp, *phsLp, *slcLp;
    LP_DIM          lp_dim;
    float           wd = 2.0; //0.5;    // smoothing width
    float           scale = 4.0;        // po2, process after this step contains FFT

    @autoreleasepool {
    printf("nciRec2\n");
        if (ac > 1) {   // for debugging
            printf("arg = %s\n", av[1]);
            asc_path = [NSString stringWithFormat:@"meas_files/meas.asc_%s", av[1]];
            out_path = [NSString stringWithFormat:@"meas_files/meas.out_%s", av[1]];
            img = [RecImage imageWithMeasAsc:asc_path andMeasData:out_path lpDim:&lp_dim];
        } else {
            img = [RecImage imageWithMeasAsc:@"meas.asc" andMeasData:@"meas.out" lpDim:&lp_dim];
        }

        system("rm IMG_*");
    // read raw / ft / phase correction / registration
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_nci2" andMeasData:@"meas_files/meas.out_nci2" lpDim:&lp_dim];  // (Siemens)
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_nci7" andMeasData:@"meas_files/meas.out_nci7" lpDim:&lp_dim];  // (Siemens)
    //    img = [RecImage imageWithMeasAsc:@"meas_files/meas.asc_0508" andMeasData:@"meas_files/meas.out_0508" lpDim:&lp_dim];  // (Siemens)

 //       [img removePointLoops];
        // kx : nSamples    
        // ky : nLin
        // kz : nSlc
        // avg: nAcq
        // phs: nRep
        avgLp = [RecLoop findLoop:@"avg"];
        phsLp = [RecLoop findLoop:@"phs"];
        slcLp = [RecLoop findLoop:@"kz"];

        [img fft2d:REC_FORWARD];
        [img freqCrop];
        [img epiPcorr];
//[img saveAsKOImage:@"IMG_epi_corr"];      // even-odd correction
        [img baseline];     // (NCI)
//[img saveAsKOImage:@"IMG_porr_fine"];       // fine phase correction
        [img rotImage:lp_dim.inplaneRot];   // Siemens
        [img reorderSlice]; // (Siemens)

    //    img = [img rotShiftForEachOfLoop:slcLp scale:scale]; // wobble
        img = [img scaleXYBy:scale];    // simple oversampling

        [img saveAsKOImage:@"IMG_base"];    // straight recon [phs avg kz ky kx]
 
        avg = [img avgForLoop:avgLp];       // avg for "avg" loop
        avg = [avg avgForLoop:phsLp];
//        [avg saveAsKOImage:@"IMG_avg"];    // straight recon [phs avg kz ky kx]

    // remove real part
        [img removeRealPart];
    // k-filter
    //    ref = [img copy];
    //    [ref changeLoop:phsLp dataLength:1 offset:0];  // first phs
     //   [ref removePointLoops];                         // remove ph loop

        ref = [img sliceAtIndex:0 forLoop:phsLp];
        ref = [ref kFilt:1];      // ## fix direction (use def of paper version)
        [ref gauss2DLP:wd]; // smoothing

   //     [ref saveAsKOImage:@"IMG_ref"];
   //     phs = [img copy];
   //     [phs changeLoop:phsLp dataLength:1 offset:1];   // take 2nd phase
    //    [phs changeLoop:phsLp dataLength:1 offset:2];   // take 3rd phase
    //    [phs subImage:ref];

    //    phs = [img removeSliceAtIndex:0 forLoop:phsLp];
		phs = [img copy];
		[phs removeSliceAtIndex:0 forLoop:phsLp];

        phs = [phs kFilt:1];      // ## fix direction (use def of paper version)
        [phs gauss2DLP:wd]; // smoothing
	
    //    [phs saveAsKOImage:@"IMG_stm"];
    // p-image (t-test)
//        pimg = [phs pImageWithRef:ref forLoop:avgLp mode:4 ttest:YES];
    //    [pimg saveAsKOImage:@"IMG_pimg"];
    // fusion
    //    col = [avg fusionWith:pimg gain:0.02 mode:1];
    //    col = [avg fusionWith:pimg gain:0.03 mode:0];  // 0:sum, 1:over
        col = [avg fusionWith:pimg gain:0.03 mode:0x00];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
        [col saveAsKOImage:@"IMG_fus1"];
        col = [col oversample];
        [col saveAsKOImage:@"IMG_fus1_hi"];
        col = [avg fusionWith:pimg gain:0.03 mode:0x01];  // 00:sum, 10:over, 00:p-mag, 01:p-dir
        [col saveAsKOImage:@"IMG_fus2"];
        col = [col oversample];
        [col saveAsKOImage:@"IMG_fus2_hi"];

        avg = [avg oversample];
        [avg saveAsKOImage:@"IMG_avg_hi"];
        pimg = [pimg oversample];
        [pimg saveAsKOImage:@"IMG_pimg_hi"];

    }   // @autoreleasepool
    return 0;
}

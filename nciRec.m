//
//	NCI basic recon
//	-> reproduce functionarity of C version nci_rec.c
//  not used any more (nciRec2 reads meas.out directly)
//

#import <RecKit/RecKit.h>

int
main()
{
	RecImage        *img, *slice;
    RecLoop         *slcLp;
    NSMutableArray  *loops;
    RecLoopControl  *lc;
    RecLoopIndex    *li;
    LP_DIM          lp_dim;
    int             i, n;
    NSString        *path;

    printf("nciRec\n");

    @autoreleasepool {
    // read raw
        img = [RecImage imageWithMeasAsc:@"meas.asc" andMeasData:@"meas.out" lpDim:&lp_dim];
        [img fft2d:REC_INVERSE];
        [img freqCrop];
        [img epiPcorr];
    //[img saveAsKOImage:@"IMG_tmp"];
        slcLp = [RecLoop findLoop:@"kz"];
        loops = [NSMutableArray arrayWithArray:[img loops]];
        [loops removeObject:slcLp];
        slice = [RecImage imageOfType:[img type] withLoopArray:loops];

    //  FFT to make corrected raw
        [img fft2d:REC_FORWARD];
        n = [slcLp dataLength];
        if (n == 0) {
            printf("Slice loop not found\n");
            exit(0);
        }
    //printf("n_slice = %d\n", n);
        lc = [img control];
        li = [lc loopIndexForLoop:slcLp];
        [lc rewind];
        [lc deactivateLoop:slcLp];
        for (i = 0; i < n; i++) {
            printf("slice %d\n", i);
            [li setCurrent:i];
            [slice copyImage:img withControl:lc];
            path = [NSString stringWithFormat:@"raw%03d.block", i];
            [slice saveAsKOImage:path];
        }
    }   // autoreleasepool
	return 0;
}

//
//  "wobbling"
//
// === plans ===
//  - single shift / rot
//  - avg of multiple images with different shift / rot
//

#import <RecKit/RecKit.h>

RecImage    *wobble(RecImage *phantom, int n, float sft);   // freq domain shift
RecImage    *wobble2(RecImage *phantom, int n);  // space domain rot / shift

int
main()
{
    RecImage    *ph;    // phantom
    RecImage    *img;   // sample
    RecImage    *ft;
    int         nAvg = 20;

    ph = [RecImage imageWithKOImage:@"phantom/PH_res.img"]; // 1024 x 1024

//    img = wobble(ph, nAvg, 10.0);
    img = wobble2(ph, nAvg);

    [img saveAsKOImage:@"phantom/PH_wobble.img"];
    img = [img rotShiftForEachOfLoop:nil scale:2.0];
    [img saveAsKOImage:@"phantom/PH_wobble_sft.img"];
    ft = [img copy];
    [ft fft2d:REC_INVERSE];
    [ft saveAsKOImage:@"phantom/PH_sft_ft.img"];
    
    img = [img avgForLoop:[img zLoop]];
    [img saveAsKOImage:@"phantom/PH_wobble_avg.img"];
    ft = [img copy];
    [ft fft2d:REC_INVERSE];
    [ft saveAsKOImage:@"phantom/PH_avg_ft.img"];

    return 0;
}

RecImage *
wobble(RecImage *ph, int n, float maxShift)
{
    RecImage        *img;   // result
    RecImage        *sft;   // tmp image for shifting
    RecLoop         *avgLp;   // avg loop
    RecLoopControl  *srcLc, *dstLc;
    RecLoopIndex    *li;
    int             i, j, k;
    int             xDim = [ph xDim];
    int             yDim = [ph yDim];
    float           dx, th, cs, sn;
    float           *p, *q, re, im;

    [ph fft2d:REC_INVERSE];
//[ph saveAsKOImage:@"phantom/PH_res_k.img"];

// make image array
    img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:256 yDim:256 zDim:n];
    [img setUnit:REC_FREQ forLoop:[img xLoop]];
    [img setUnit:REC_FREQ forLoop:[img yLoop]];
    avgLp = [img zLoop];
    dstLc = [img control];
    [dstLc rewind];
    [dstLc deactivateLoop:avgLp];
    li = [dstLc loopIndexForLoop:avgLp];

// fft / shift / ift
    for (k = 0; k < n; k++) {
        [li setCurrent:k];
        sft = [ph copy];
        p = [sft data];
        q = p + [sft dataLength];
        // ft base shift
        dx = -(float)k * maxShift / n;
        for (i = 0; i < yDim; i++) {
            for (j = 0; j < xDim; j++) {
                th = ((float)j - xDim/2) / xDim * dx * M_PI * 2.0;
				cs = cos(th);
				sn = sin(th);
				re = p[i * xDim + j];
				im = q[i * xDim + j];
				p[i * xDim + j] = re * cs + im * sn;
				q[i * xDim + j] =-re * sn + im * cs;
            }
        }

        [sft crop:[sft xLoop] to:256];
        [sft crop:[sft yLoop] to:256];
        srcLc = [sft control];
        [img copyImage:sft dstControl:dstLc srcControl:srcLc];
    }
    [img fft2d:REC_FORWARD];
    return img;
}

RecImage *
wobble2(RecImage *ph, int n)  // space domain rot / shift
{
    RecImage        *img;   // result
    RecImage        *param; // rot/shift amount
    RecLoop         *avgLp;   // avg loop
    int             i;
    int             xDim = [ph xDim];
    int             yDim = [ph yDim];
    float           mx;
    float           *p, *q;

// make image array
    img = [RecImage imageOfType:RECIMAGE_COMPLEX withImage:ph];
    avgLp = [RecLoop loopWithDataLength:n];
    [img replaceLoop:[img zLoop] withLoop:avgLp];
// map
	param = [RecImage imageOfType:RECIMAGE_MAP withLoops:avgLp, nil];	// array of 2d vectors
	p = [param data];
	q = p + [param dataLength];

    [img copyImage:ph];
//[img saveAsKOImage:@"phantom/test_img"];

    // rot / shift (space)
// rotation (1D param)
    if (1) {
        mx = 0.02;
        printf("0 to %4.2f PI, %d steps\n", mx, n);
        for (i = 0; i < n; i++) {
            p[i] = M_PI * mx * i / n;
        }
        img = [img rotBy:param];
        [img saveAsKOImage:@"phantom/test_img"];
    }

// translation (2D param)
    if (1) {
        mx = 10.0;	// pixels
        printf("0 to %4.2f pixels, %d steps\n", mx, n);
        xDim = [img xDim];
        yDim = [img yDim];
        for (i = 0; i < n; i++) {
            p[i] = mx * (float)i / n / xDim;	// FOV frac (0 - 0.039)
//printf("shift = %f (pixels), %f frac\n", p[i] * xDim, p[i]);
            q[i] = p[i] * 0.3;
        }
        img = [img ftShiftBy:param];
        [img saveAsKOImage:@"phantom/test_rot"];
    }

    [img fft2d:REC_INVERSE];
    [img crop:[img xLoop] to:256];
    [img crop:[img yLoop] to:256];

    [img fft2d:REC_FORWARD];
    return img;
}


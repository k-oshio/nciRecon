//
//	combine slices
//
//  input: IMG_ica.%d, ICA_a%d.dat, IMG_mag%d
//  output: IMG_ica, IMG_mag
//
//	### not done yet
//

#import <RecKit/RecKit.h>

float       calc_corr(float *p);
RecImage    *findMaxCorr(RecImage *ss_ica, int ix);
void        rot45(RecImage *f);

int     nx, ny;

int
main(int ac, char *av[])
{
    int             i, nSl = 5;
    RecImage        *ica, *mag, *ss_ica, *ss_mag, *ica_comb;
    RecLoop         *xLp, *yLp, *slLp, *phLp;
    RecLoopControl  *lc;
    RecLoopIndex    *slLi;
    NSString        *path;

    ss_ica = [RecImage imageWithKOImage:@"IMG_ica0"];
    if (ss_ica == nil) {
        printf("image not found\n");
        exit(0);
    }
    xLp = [ss_ica xLoop];
    yLp = [ss_ica yLoop];
    phLp = [ss_ica topLoop];
    slLp = [RecLoop loopWithName:@"Slice" dataLength:5];

    ica = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:slLp, yLp, xLp, nil];
    mag = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:slLp, yLp, xLp, nil];

    ss_ica = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:phLp, yLp, xLp, nil];
    ss_mag = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:yLp, xLp, nil];

    lc = [RecLoopControl controlWithLoops:slLp, phLp, yLp, xLp, nil];

    slLi = [lc loopIndexForLoop:slLp];
    [slLi deactivate];

    for (i = 0; i < nSl; i++) {
        [slLi setCurrent:i];
    // read mag image, set to result
        path = [NSString stringWithFormat:@"IMG_mag%d", i];
        [ss_mag initWithKOImage:path];
        [mag copyImage:ss_mag withControl:lc];
    // read a matrix, calc corr, find max corr component
        path = [NSString stringWithFormat:@"ICA_a%d.dat", i];
        ica_comb = findMaxCorr(ss_ica, i);
        [ica copyImage:ica_comb withControl:lc];
    }
    [ica saveAsKOImage:@"IMG_ica"];
    [mag saveAsKOImage:@"IMG_mag"];
// make 45deg rotated version (do this in kFilt)
    rot45(ica);
    [ica saveAsKOImage:@"IMG_ica_ob"];
}

RecImage *
findMaxCorr(RecImage *ss_in, int ix)
{
    RecImage        *comb_img, *tmp_img;
    RecLoop         *xLp, *yLp;
    RecLoopControl  *lc;
    RecLoopIndex    *li;
    NSString        *path, *content, *str;
    NSArray         *lines;
    NSError         *err;
    NSScanner       *scanner;
    int             i, j;
    float           *buf, *corr;
    float           thres = 0.30;//0.40;//30;

    path = [NSString stringWithFormat:@"ICA_a%d.dat", ix];
    content = [NSString stringWithContentsOfFile:path
                encoding:NSASCIIStringEncoding error:&err];
    lines = [content componentsSeparatedByString:@"\n"];
    nx = [lines count] - 1;
    str = [lines objectAtIndex:0];
    scanner = [NSScanner scannerWithString:str];
    for (i = 0; ; i++) {
        if (![scanner scanFloat:nil]) break;
    }
    ny = i - 1;

    buf = (float *)malloc(ny * nx * sizeof(float));
    corr = (float *)malloc(ny * sizeof(float));

    for (i = 0; i < nx; i++) {
        str = [lines objectAtIndex:i];
        scanner = [NSScanner scannerWithString:str];
        [scanner scanFloat:nil];    // index
        for (j = 0; j < ny; j++) {
			[scanner scanFloat:buf + j * nx + i];
        }
    }

// calc corr
    for (i = 0; i < ny; i++) {
        corr[i] = calc_corr(buf + i * nx);
    }
// make comb_img
    path = [NSString stringWithFormat:@"IMG_ica%d", ix];
    [ss_in initWithKOImage:path];
    xLp = [ss_in xLoop];
    yLp = [ss_in yLoop];
    comb_img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:yLp, xLp, nil];
    tmp_img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:yLp, xLp, nil];
    lc = [ss_in control];
    li = [lc topLoopIndex];
    [li deactivate];

    [comb_img clear];
    for (i = 0; i < ny; i++) {
        if (fabs(corr[i]) > thres) {
            [li setCurrent:i];
            [tmp_img copyImage:ss_in withControl:lc];
            [tmp_img multByConst:corr[i]];
            [comb_img addImage:tmp_img];
        }
    }
    free(buf);
    free(corr);

    return comb_img;
}

// extract and delete
int
plotCsv(char *csvPath, char *datPath)
{
    FILE    *fp;
	int		i, j, ix;
    int     xdim, ydim;
    int     c;
    float   *buf, *p, val;
    float   len, ip, mn, mx;

    fp = fopen(csvPath, "r");

    xdim = 0;
    while (((c = getc(fp)) != EOF) && (c != '\n')) {
        if (c == ',') xdim++;
    }
    xdim += 1;

    rewind(fp);
    ydim = 0;
    while ((c = getc(fp)) != EOF) {
        if (c == '\n') ydim++;
    }
    
    printf("xdim = %d, ydim = %d\n", xdim, ydim);

// read csv
    buf = (float *)malloc(xdim * ydim * sizeof(float));
    rewind(fp);
    p = buf;
    for (i = 0; i < ydim; i++) {
        for (j = 0; j < xdim; j++) {
            fscanf(fp, "%f,", p + j);
         }
         fscanf(fp, "\n");
         p += xdim;
    }
    fclose(fp);

// write .dat
    fp = fopen(datPath, "w");
    for (i = 0; i < ydim; i++) {
        fprintf(fp, "%d ", i);
        for (j = 0; j < xdim; j++) {
            val = *(buf + i * xdim + j);
            fprintf(fp, "%7.4f ", val);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

// find max corr with stim
    for (j = 0; j < xdim; j++) {
        mn = 0;
        for (i = 0; i < ydim; i++) {
            val = *(buf + i * xdim + j);
            mn += val;
        }
        mn /= ydim;
        for (i = 0; i < ydim; i++) {
            *(buf + i * xdim + j) -= mn;
        }
    }
    mx = 0;
    for (j = 0; j < xdim; j++) {
        ip = len = 0;
        for (i = 0; i < ydim; i++) {
            val = *(buf + i * xdim + j);
            if (i % 2 == 0) {
                ip += val;
            } else {
                ip -= val;
            }
            len += val*val;
        }
        ip /= sqrt(len * ydim);
        printf("%d %f\n", j, ip);
        if (fabs(mx) < fabs(ip)) {
            mx = ip;
            ix = j;
        }
    }
    printf("max is %f at %d\n", mx, ix);
    free(buf);

    return ix;
}

float
calc_corr(float *p)
{
    float   c, len, val;
    int     i;
    
    c = len = 0;
    for (i = 0; i < nx; i++) {
        val = p[i];
        if (i % 2 == 0) {
            c += val;
        } else {
            c -= val;
        }
        len += val*val;
    }
    c /= sqrt(len * nx);

    return c;
}

void
rot45(RecImage *f)
{
    float   *p, *q;
    float   re, im;
    int     i, len;

    len = [f dataLength];
    p = [f data];
    q = p + len;

    for (i = 0; i < len; i++) {
        re = p[i] + q[i];
        im = p[i] - q[i];
        p[i] = re;
        q[i] = im;
    }
}



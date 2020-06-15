//

#import <RecKit/RecKit.h>

void    findMaxCorr(NSString *path);
int     nx, ny;
float   calc_corr(float *p);

int
main()
{
    @autoreleasepool {

        findMaxCorr(@"ICA_a2.dat");
    }
    return 0;
}
/*
0 -0.706581
1 -0.398663
2 0.430383
3 0.461998
4 -0.302155
5 -0.076466
6 -0.047893
7 0.062624
8 0.207787
9 -0.150841
10 -0.385293
11 0.035257
12 -0.072136
13 -0.153658
14 -0.479178
*/
void
findMaxCorr(NSString *path)
{
    NSString    *content, *str;
    NSArray     *lines;
    NSError     *err;
    NSScanner   *scanner;
    int         i, j;
    float       *buf, *corr, *comb;
    float       thres = 0.30;

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
    comb = (float *)malloc(nx * sizeof(float));
//    printf("n = %d, nc = %d\n", n, nc);
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
// pickup ones above thres
    for (j = 0; j < nx; j++) {
        comb[j] = 0;
    }
    for (i = 0; i < ny; i++) {
    //    printf("%d %f\n", i, corr[i]);
        if (fabs(corr[i]) > thres) {
         //   printf("%d %f\n", i, corr[i]);
            for (j = 0; j < nx; j++) {
                comb[j] += buf[i * nx + j] * corr[i];
            }
        }
    }
    for (j = 0; j < nx; j++) {
        printf("%d %f\n", j, comb[j]);
    }

    free(buf);
    free(corr);
    free(comb);
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

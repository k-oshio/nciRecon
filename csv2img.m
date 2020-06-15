//
//	csv to RecImage conversion (with mask)
//  (for Octave processing)
//
//  input: IMG_s.csv, IMG_mask
//  output: IMG_ica
//
//	### not done yet
//

#import <RecKit/RecKit.h>

int
main(int ac, char *av[])
{
	RecImage	*f, *mask;
	RecLoop		*kx, *ky, *kz;
	float		*p, *q, *msk;
	int			i, j, ix, c;
    int         imgSize, nSlice, nEl, dataLength;
    FILE        *fp;
    NSArray     *lines, *items;
    NSString    *content;
    NSError     *err;

//    content = [NSString stringWithContentsOfFile:@"IMG_s.csv"
    content = [NSString stringWithContentsOfFile:@"IMG_x.csv"
        encoding:NSASCIIStringEncoding error:&err];
    lines = [content componentsSeparatedByString:@"\n"];
    nSlice = (int)[lines count];

    for (i = 0; i < nSlice; i++) {
        content = [lines objectAtIndex:i];
        items = [content componentsSeparatedByString:@","];
        printf("Line:%d, nItem:%lu\n", i, [items count]);
    }
exit(0);




    fp = fopen("IMG_s.csv", "r");
    // count # of lines
    nSlice = 0;
    while ((c = getc(fp)) != EOF) {
        if (c == '\n') nSlice++;
    }
    mask = [RecImage imageWithKOImage:@"IMG_mask"];
    msk = [mask data];
    dataLength = [mask dataLength];

    nEl = 0;
    for (i = 0; i < dataLength; i++) {
        if (msk[i] != 0) nEl++;
    }
    printf("nSlice = %d, nEl = %d\n", nSlice, nEl);

	kx = [RecLoop loopWithName:@"kx" dataLength:128];
	ky = [RecLoop loopWithName:@"ky" dataLength:128];
	kz = [RecLoop loopWithName:@"kz" dataLength:nSlice];
	f = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:kz, ky, kx, nil];
    dataLength = [f dataLength];
    imgSize = 128 * 128;

    p = [f data];
    q = p + [f dataLength];
    msk = [mask data];

    rewind(fp);
    for (i = 0; i < nSlice; i++) {
        for (j = ix = 0; j < imgSize; j++) {
            if (msk[j] != 0) {
                fscanf(fp, "%f,", &p[j]);
                ix++;
            }
            if (ix >= nEl) break;
        }
        for (j = ix = 0; j < imgSize ; j++) {
            if (msk[j] != 0) {
                fscanf(fp, "%f,", &q[j]);
                ix++;
            }
           if (ix >= nEl) break;
        }
        fscanf(fp, "0\n");
        p += imgSize;
        q = p + dataLength;
    }
    fclose(fp);

	[f saveAsKOImage:@"IMG_ica"];
}

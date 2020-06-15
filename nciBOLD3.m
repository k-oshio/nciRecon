//
//	nciBOLD2
//	3-12-2018		forked off from nciRecBOLD1
//					phase of GE EPI seq (BOLD set, or single stim set)
//	=== plans ===
//	just look at phase part of BOLD after unwrapping along time direction
//


#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

int
main(int ac, char *av[])
{
    RecImage        *img, *img_a, *img_u;
	RecImage		*img_f, *avg;
	RecImage		*mask;
	RecLoop			*timeLp;
	int				nTime;
	int				xDim, yDim;
	// PCA
	Num_mat			*A;
	Num_svd_result	*sres;
	// locations
	NSString		*base = @"/Users/oshio/images/NCI/NIRS/2018-0611/results";
	int				series = 3;
	NSString		*path;

    printf("nciRecBOLD3 (6-14-2018)\n");

    @autoreleasepool {
		img = nil;

// === read raw ===
		path = [NSString stringWithFormat:@"%@/%d/IMG_comb", base, series];	// KOImage with no loop structure
		img = [RecImage imageWithKOImage:path];

		[img crop:[img zLoop] to:1024 startAt:100];
		[img saveAsKOImage:@"IMG_crop"];

		timeLp = [img zLoop];
		nTime = [timeLp dataLength];
		xDim = [img xDim];
		yDim = [img yDim];

//		mask = [img sdMaskForLoop:[img zLoop] thres:0.1];	// rect mask
//		[mask hGauss2DLP:0.5];	// "Half" gaussian filter
//		[mask saveAsKOImage:@"IMG_mask"];

//		[img phase];
//		[img multByImage:mask];

//===== freq analysis
		img_f = [img copy];
		avg = [img_f avgForLoop:[img_f zLoop]];
		[img_f subImage:avg];
		[img_f fft1d:[img_f zLoop] direction:REC_INVERSE];
		[img_f saveAsKOImage:@"IMG_f"];

//====== PCA
		img_a = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim * yDim yDim:nTime];
		[img_a copyImageData:img];
		A = Num_im_to_m(img_a);
		sres = Num_svd(A);
		img_u = Num_m_to_im(sres->U);
		[img_u trans];
		[img_u crop:[img_u yLoop] to:32 startAt:0];
		[img_u saveAsKOImage:@"IMG_U"];
		img_a = Num_m_to_im(sres->Vt);
		[img copyImageData:img_a];
		[img saveAsKOImage:@"IMG_PCA"];

		Num_free_svd_result(sres);
		Num_free_mat(A);
		

    }   // @autoreleasepool

    return 0;
}


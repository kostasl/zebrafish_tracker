#ifndef CONFIG_H
#define CONFIG_H


#define ZTF_FISHCONTOURSIZE          45
#define ZTF_TAILFITMAXITERATIONS     10 //For Spine To Contour Tail Fitting
#define ZTF_TAILSPINECOUNT          8
#define USE_CUDA
#undef  USE_CUDA_FOR_TEMPLATE_MATCHING //Can Be Very Slow When Using Current Algorithm for searching through Cache
#define TEMPLATE_COL_SEARCH_REGION  15 //Optimization How Far to Search Along Col in the Cache From the Starting POint

#endif // CONFIG_H

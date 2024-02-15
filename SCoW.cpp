
/***

Reference : Watershed Superpixels, Zhongwen Hu, Qin Zou, Qingquan Li, IEEE ICIP 2015

The code is an extention of the watershed superpixel described in the reference.

****/

#include <mex.h>

#include <vector>
#include <queue>
#include <math.h>
#include <set>
#include <map>
#include <string>
#include "stdio.h"

typedef struct {
	int x;
	int y;
}zwPoint;

const int INIT = -3;
const int MASK = -2;
const int WATERSHED = -1;

#define Num_of_Queues  750

double factorial(double num)
{
	if (num==0 || num==1)
		return 1;
	   return (num * factorial(num - 1));
}
//ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ë¹ï¿½ï¿½ï¿½ï¿½Ë²ï¿½ï¿½ï¿?
void create_filters(int WinSize, double *smofil_, double *diffil_)
{
	int i;
	double w;
	int WL_ = WinSize/2;
	for (i=-WL_; i<=WL_; i++){
		w = exp(double(-2*WL_))*factorial(2*WL_)/(factorial(WL_-i)*factorial(WL_+i) );
		smofil_[i+WL_] = w;
		diffil_[i+WL_] = (2*i*w)/WL_;
	}
	double t1 = 0.0f;
	double t2 = 0.0f;
	for (i=0;i<WinSize;i++){
		t1 = smofil_[i]>t1?smofil_[i]:t1;
		t2 = diffil_[i]>t2?diffil_[i]:t2;
	}
	for (i=0;i<WinSize;i++){
		smofil_[i] /= t1;
		diffil_[i] /= t2;
	}
	
}

void gaussian_diff_filter(double *pImage, float *pGrad, int ImgSizeX, int ImgSizeY, int WinSize)
{
	
	double* sf = new double[WinSize]; //smooth filter
	double* df = new double[WinSize]; //diff filter
	//  unsigned char* im;
	
	float* tim;
	double sum = 0;
	double sum1 = 0;
	int i, j, k;
	
	//create kernels
	create_filters(WinSize,sf,df);
	int m_nWidth = ImgSizeX;
	int m_nHeight = ImgSizeY;
	
	float *pGradX = new float[m_nWidth*m_nHeight];
	float *pGradY = new float[m_nWidth*m_nHeight];
	
	//  im = cim->im_;
	tim = new float[m_nWidth*m_nHeight];
	for (i=0; i<m_nWidth*m_nHeight; i++){
		pGradX[i] = pGradY[i] = 0;
		tim[i] = pImage[i];
	}
	int WL_ = WinSize/2;
	//filter image x
	//smooth on y
	for (i=0; i<m_nWidth; i++){
		for (j=WL_; j<(m_nHeight-WL_); j++){
			sum = 0;
			for (k=-WL_; k<=WL_; k++)
				sum += (float)sf[k+WL_]*pImage[(j+k)*m_nWidth+i];
			tim[j*m_nWidth+i] = sum;
		}
	}
	//diff on x
	for (j=0; j<m_nHeight; j++){
		for (i=WL_; i<(m_nWidth-WL_); i++){
			sum = 0;
			for (k=-WL_; k<=WL_; k++)
				sum += df[k+WL_]*tim[j*m_nWidth+i+k];
			pGradX[j*m_nWidth+i] = (float) (sum);
		}
	}
	
	//filter image y
	for (i=0; i<m_nWidth*m_nHeight;i++)
		tim[i] = pImage[i];
	//  im = cim->im_;
	//smooth on x
	for (j=0; j<m_nHeight; j++){
		for (i=WL_; i<(m_nWidth-WL_); i++){
			sum = 0;
			for (k=-WL_; k<=WL_; k++)
				sum += (float)sf[k+WL_]*pImage[j*m_nWidth+i+k];
			tim[j*m_nWidth+i] = sum;
		}
	}
	
	//diff on y
	for (i=0; i<m_nWidth; i++){
		for (j=WL_; j<(m_nHeight-WL_); j++){
			sum = 0;
			for (k=-WL_; k<=WL_; k++)
				sum += df[k+WL_]*tim[(j+k)*m_nWidth+i];
			pGradY[j*m_nWidth+i] = (float) (sum);
		}
	}
	double weight = 1.0/(0.3333333 * WinSize);
	double tmp;
	for (i=0;i<m_nWidth*m_nHeight;i++){
		tmp = pow((double)(pGradX[i]*pGradX[i] + pGradY[i]*pGradY[i]),0.5);
		pGrad[i] = weight*tmp<255? weight*tmp: 255;	}
	delete []pGradX;
	delete []pGradY;
	delete [] tim;
	delete []sf;
	delete []df;
}
void get_gaussian_gradient(double *R, double *G, double *B, int iWidth, int iHeight, int WinSize,float *deltar)
{
	float *deltar_r = new float[iWidth * iHeight];
	float *deltar_g = new float[iWidth * iHeight];
	float *deltar_b = new float[iWidth * iHeight];

	gaussian_diff_filter(R,deltar_r,iWidth,iHeight,WinSize);
 	gaussian_diff_filter(G,deltar_g,iWidth,iHeight,WinSize);
 	gaussian_diff_filter(B,deltar_b,iWidth,iHeight,WinSize);
	int i;
	for(i=0;i<iWidth * iHeight;i++){
		deltar[i] = (deltar_r[i]+deltar_g[i]+deltar_b[i]) * 0.33333f * 255;	}
	delete []deltar_r; deltar_r = NULL;
	delete []deltar_g; deltar_g = NULL;
	delete []deltar_b; deltar_b = NULL;
}
void set_random_seeds(int *label, float *deltar, int iwidth, int iheight, int nSeeds, std::vector<zwPoint> &seedPos)
{
    int i;
    int SuperNum = 0;
    zwPoint tmpFPoint;
    for(i=0;i<nSeeds;i++)
    {
        int CentroX = 2+rand()%(iwidth-4);
        int CentroY = 2+rand()%(iheight-4);
		label[CentroY*iwidth+CentroX] = SuperNum;\
		tmpFPoint.y = CentroY;
		tmpFPoint.x = CentroX;
		seedPos.push_back(tmpFPoint); 
		SuperNum = SuperNum + 1;
    }

}
void set_uniform_seeds(int *label, float *deltar, int iwidth, int iheight, int winsize, std::vector< zwPoint > &seedPos)
{
	int i;
    int tmppos;
	int seedpos;
	int _WinSize = winsize;
	int halfSize = _WinSize/2;
    int imagelen = iwidth*iheight;
	if (winsize<7){
		winsize = 7;}
	for ( i=0; i<imagelen; i++) {
		label[i] = -3;  }
	int x,y;
	float minGradient = 1000;
	
	int blockSizeX = iwidth/_WinSize;
	float xDist = ((float)iwidth)/blockSizeX;
	int blockSizeY = iheight/_WinSize ;
	float yDist = ((float)iheight)/blockSizeY;
	
	float halfStepX = xDist*0.5f;
	float halfStepY = yDist*0.5f;
	
	int k=0;
	int j;
	
	int CentroX,CentroY;
	zwPoint tmpFPoint;
	int SuperNum = 0;
    //block
	for (j = 0; j < blockSizeY; j++){
		for (i = 0; i < blockSizeX; i++){
			CentroY = (0.5+j)*yDist + rand()%3-1;
			CentroX = (0.5+i)*xDist + rand()%3-1;
			minGradient = 1000;
			for (y=CentroY-2;y<CentroY+2;y++){
				for (x=CentroX-2;x<CentroX+2;x++){
					tmppos = y*iwidth+x;
					if (deltar[tmppos] < minGradient){
						minGradient = deltar[tmppos];
						seedpos = tmppos;
					}
				}
			}
			label[CentroY*iwidth+CentroX] = SuperNum;
			tmpFPoint.y = seedpos/iwidth;
			tmpFPoint.x = seedpos-iwidth*tmpFPoint.y;
			seedPos.push_back(tmpFPoint); 
			SuperNum = SuperNum + 1;
		}
	}
    /*
    		for (j = 0; j < blockSizeY; j++)
		{
			int off1 = j%2;
			float off2 = off1 * halfStepX ;
			for (i = 0; i < blockSizeX-off1; i++)
			{
				CentroY = (0.5+j)*yDist + rand()%3 -1;
				CentroX = (0.5+i)*xDist + off2 + rand()%3 - 1;
				tmppos = CentroY*xDist + CentroX;
				minGradient = 500;
				for (y=CentroY-2;y<CentroY+2;y++)
				{
					for (x=CentroX-2;x<CentroX+2;x++)
					{
						tmppos = y*iwidth+x;
						if (deltar[tmppos] < minGradient)
						{
							minGradient = deltar[tmppos];
							seedpos = tmppos;
						}
					}
				}
				label[seedpos] = SuperNum;
				tmpFPoint.y = seedpos/iwidth;
				tmpFPoint.x = seedpos-iwidth*tmpFPoint.y;
				seedPos.push_back(tmpFPoint); //¼ÇÂ¼ÖÖ×Ó±ê¼ÇµÄÎ»ÖÃ
				SuperNum++;
			}
		}*/

}
void ExpandBasin2(int center,int pos, int *label)
{
    
    if ( label[pos]>WATERSHED )
    {
		if (label[center]<WATERSHED){
            label[center] = label[pos];}
		else if (label[center]==WATERSHED){	
            }
		else if (label[center]!=label[pos]){
			label[center] = WATERSHED; }
    }
    
}
void superpixel_flood(double lambda,/*double beta, */std::vector< zwPoint > m_SeedPos, float *deltar,int Winsize, int ImgSizeX, int ImgSizeY, int *label)
{
	int i;
    int Imglen = ImgSizeX * ImgSizeY;
    for (i=0;i<ImgSizeX;i++)
    {
		label[i] = -1;
		label[Imglen-1-i] = -1;
    }
    for (i=0;i<ImgSizeY;i++)
    {
		label[i*ImgSizeX] = -1;
		label[(i+1)*ImgSizeX-1] = -1;
    }
	double  mean = 0;
	for(i=0;i<Imglen;i++)
	{
		mean += deltar[i];
	}
	mean = mean/Imglen;

	double beta;

	beta = pow(mean, (double)0.75f);

	mexPrintf("Automatical Estimated Beta = %.3f\nlambda = %.3f\n",beta,lambda);
	
	float m_bShapeConstraint[Num_of_Queues];
    for (i=0;i<Num_of_Queues;i++)
        m_bShapeConstraint[i] = 0;

	int ObjNum = m_SeedPos.size();
	mexPrintf("Num.of.Superpixels %d\n", ObjNum);

	double maxValue = 0;

	for (i=0;i<Num_of_Queues;i++)
	{
		m_bShapeConstraint[i] = exp(-(double)i/beta)*lambda;
	}
    double *distweight = new double[500];
    for(i=0;i<500;i++)
    {
        double tmp = -sqrt(i);
        distweight[i] = 1-exp(tmp*0.25);
    }
	unsigned char *isProcessed = new unsigned char[ImgSizeX*ImgSizeY];
    //memset(isProcessed,0,ImgSizeX*ImgSizeY*sizeof(unsigned char));
    for(i=0;i<ImgSizeX*ImgSizeY;i++)
        isProcessed[i] = 0;
    std::queue< zwPoint > PriQue[Num_of_Queues];//
    int j;
    int TmpY,TmpPos;
    zwPoint TmpPt;
    for (i=1;i<ImgSizeY-1;i++)
    {
		TmpY = i*ImgSizeX;
		for (j=1;j<ImgSizeX-1;j++)
		{
			TmpPos = TmpY + j;
			if (label[TmpPos]<WATERSHED)
				if(label[TmpPos-1]>WATERSHED
					||label[TmpPos+1]>WATERSHED
					||label[TmpPos+ImgSizeX]>WATERSHED
					||label[TmpPos-ImgSizeX]>WATERSHED
					)
				{
					TmpPt.x = j; TmpPt.y = i;
                    PriQue[0].push(TmpPt);                    
					isProcessed[TmpPos] = 1;
				}
		}
    }
    zwPoint TmpPt2;
    int time = 0;
    unsigned short Index[Num_of_Queues]; 
    for (i=0;i<Num_of_Queues;i++)
    {
		Index[i] = i;
    }
    for (i=0;i<Num_of_Queues;i++)
    {   
		
		int k=0;
		for (k=0;k<i;k++)
		{
			Index[k] = i;
		}
		
		double dis2seed = 0.0f;
		zwPoint seedPos;
		while(!PriQue[i].empty())
		{
			TmpPt = PriQue[i].front();
			PriQue[i].pop();
			TmpPos = TmpPt.y*ImgSizeX + TmpPt.x;
			ExpandBasin2(TmpPos,TmpPos-1,label);        
			ExpandBasin2(TmpPos,TmpPos+1,label);        
			ExpandBasin2(TmpPos,TmpPos-ImgSizeX,label); 
			ExpandBasin2(TmpPos,TmpPos+ImgSizeX,label); 
			if (label[TmpPos]>=0)
			{
				seedPos.y = m_SeedPos[label[TmpPos]].y;
				seedPos.x = m_SeedPos[label[TmpPos]].x;

				if (label[TmpPos-1]<WATERSHED&&!isProcessed[TmpPos-1])
				{
					TmpPt2.x = TmpPt.x-1;TmpPt2.y = TmpPt.y;
					dis2seed = (seedPos.x-TmpPt2.x)*(seedPos.x-TmpPt2.x)+(seedPos.y-TmpPt2.y)*(seedPos.y-TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
                     int p = Index[int(deltar[TmpPos-1]*distweight[int(dis2seed)])] + m_bShapeConstraint[i]*dis2seed ;
					PriQue[p].push(TmpPt2);
					isProcessed[TmpPos-1] = 1;
				}
				if (label[TmpPos+1]<WATERSHED&&!isProcessed[TmpPos+1])//ï¿½ï¿½
				{
					TmpPt2.x = TmpPt.x+1;TmpPt2.y = TmpPt.y;
					dis2seed = (seedPos.x-TmpPt2.x)*(seedPos.x-TmpPt2.x)+(seedPos.y-TmpPt2.y)*(seedPos.y-TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
                   int p = Index[int(deltar[TmpPos+1]*distweight[int(dis2seed)])] + m_bShapeConstraint[i]*dis2seed ;
					PriQue[p].push(TmpPt2);
					isProcessed[TmpPos+1] = 1;
				}
				if (label[TmpPos-ImgSizeX]<WATERSHED&&!isProcessed[TmpPos-ImgSizeX])//ï¿½ï¿½
				{
					TmpPt2.x = TmpPt.x;TmpPt2.y = TmpPt.y-1;
					dis2seed = (seedPos.x-TmpPt2.x)*(seedPos.x-TmpPt2.x)+(seedPos.y-TmpPt2.y)*(seedPos.y-TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
                    int p = Index[int(deltar[TmpPos-ImgSizeX]*distweight[int(dis2seed)])] + m_bShapeConstraint[i]*dis2seed ;
					PriQue[p].push(TmpPt2);
					isProcessed[TmpPos-ImgSizeX] = 1;
				}
				if (label[TmpPos+ImgSizeX]<WATERSHED&&!isProcessed[TmpPos+ImgSizeX])//ï¿½ï¿½
				{
					TmpPt2.x = TmpPt.x;TmpPt2.y = TmpPt.y+1;
					dis2seed = (seedPos.x-TmpPt2.x)*(seedPos.x-TmpPt2.x)+(seedPos.y-TmpPt2.y)*(seedPos.y-TmpPt2.y);
					dis2seed = sqrt((double)dis2seed);
                    int p = Index[int(deltar[TmpPos+ImgSizeX]*distweight[int(dis2seed)])] + m_bShapeConstraint[i]*dis2seed ;
					PriQue[p].push(TmpPt2);
					isProcessed[TmpPos+ImgSizeX] = 1;
				}
			}
		}
    }

    delete []isProcessed;
    isProcessed = NULL; 
    delete []distweight;
    distweight = NULL;
}
void findedges(int *label, unsigned char *edge, int iwidth, int iheight)
{
	int i;
	for (i=0;i<iheight*iwidth;i++){
		if (label[i]<0){
			edge[i] = 255;
		}
		else edge[i] = 0;
	}
}
double distance(double *R, double *G, double *B, int pos1, int pos2)
{
	double dist = 0.0f;
	dist += (R[pos1]-R[pos2])*(R[pos1]-R[pos2]);
	dist += (G[pos1]-R[pos2])*(G[pos1]-R[pos2]);
	dist += (B[pos1]-R[pos2])*(B[pos1]-R[pos2]);
	return dist;
}
void ProcEdgePixs(int *label, double *R, double *G, double *B, std::vector< zwPoint > seedPos,int ImgSizeX, int ImgSizeY)
{
    int x[4]={-1,1,-ImgSizeX,ImgSizeX};
    int xx[4]={-1,1,0 ,0};
    int yy[4]={0 ,0,-1,1};
    int i,j;
    int xstart,tmppos;

	unsigned char *edge = new unsigned char[ImgSizeX*ImgSizeY];
	findedges(label,edge,ImgSizeX,ImgSizeY);
    int k;
    double Dist;
    double min = 10000;
    int local;
    for (i=1;i<ImgSizeY-1;i++) {
		xstart = i*ImgSizeX;
		for (j=1;j<ImgSizeX-1;j++){
			tmppos = xstart + j;
			if (edge[tmppos]==255){
				min = 999999;
				local = 0;
				for (k=0;k<4;k++){
					if (edge[tmppos+x[k]]!=255)	{
						Dist = distance(R,G,B,tmppos,tmppos+x[k]);
						Dist += sqrt((double)((seedPos[label[tmppos+x[k]]].x-j)*(seedPos[label[tmppos+x[k]]].x-j)+(seedPos[label[tmppos+x[k]]].y-i)*(seedPos[label[tmppos+x[k]]].y-i)));
						if (Dist<min){
							local = k;
							min = Dist;
						}
					}
				}
				label[tmppos] = label[tmppos+x[local]];
				if (label[tmppos]<0){
					for (k=0;k<4;k++){
						if (label[tmppos+x[k]]>=0){
							label[tmppos] = label[tmppos+x[k]];
							break;
						}
					}
				}
				
			}
		}
    }
    for (i=0;i<ImgSizeY;i++) {
		label[i*ImgSizeX] = label[i*ImgSizeX+1];
		label[i*ImgSizeX+ImgSizeX-1] = label[i*ImgSizeX+ImgSizeX-2];
    }
    int ImgSize = ImgSizeX*ImgSizeY-ImgSizeX;
    for (i=0;i<ImgSizeX;i++) {
		label[i] = label[i+ImgSizeX];
		label[i+ImgSize] = label[i+ImgSize-ImgSizeX];
    }
	delete []edge; edge = NULL;
}

void SmootBoundaries(int *label, int ImgSizeX, int ImgSizeY)
{

    int up,down,left,right;
    int flag = 0;
    int m_b[4];
    int i,j,xstart,tmppos;
    for (i=1;i<ImgSizeY-1;i++)  {
		xstart = i*ImgSizeX;
		for (j=1;j<ImgSizeX-1;j++){
			tmppos = xstart + j;
			up = tmppos - ImgSizeX;
			down = tmppos +ImgSizeX;
			left = tmppos - 1;
			right = tmppos + 1;
			flag = 0;
			m_b[0] = m_b[1] = m_b[2] = m_b[3] = 0;
			
			if (label[tmppos]==label[up]) { flag++; m_b[0] = 1;}
			if (label[tmppos]==label[left]) { flag++; m_b[1] = 1;}
			if (label[tmppos]==label[down]) { flag++; m_b[2] = 1;}
			if (label[tmppos]==label[right]) { flag++; m_b[3] = 1;}
			if (flag==1){
				if (label[up]==label[down]||label[left]==label[right])				{
					if(!m_b[0]) { label[tmppos] = label[up]; }
					else if(!m_b[1]) { label[tmppos] = label[left]; }
					else if(!m_b[2]) { label[tmppos] = label[down]; }
					else if(!m_b[3]) { label[tmppos] = label[right]; }
				}
			}
		}
    }
}
void SCoW_rgb( double *R, double *G, double *B, int *label, unsigned char *boundary, int iwidth, int iheight, double winSize, double lambda)
{
	
	float *deltar = new float[iwidth*iheight];

	int imgsize = iwidth * iheight;
	
	get_gaussian_gradient(R,G,B,iwidth,iheight,5,deltar);

	mexPrintf("Compute gradient!\n");
    
    std::vector< zwPoint > seedPos;

 	set_uniform_seeds(label,deltar,iwidth, iheight, winSize, seedPos);
    //set_random_seeds(label,deltar,iwidth,iheight,1000,seedPos);

	mexPrintf("Set Uniform Seeds!\n");

 	superpixel_flood(lambda,seedPos,deltar,winSize,iwidth,iheight,label);
 
	mexPrintf("Flood!\n");


	delete []deltar; deltar = NULL;
    
    for(int m=0;m<iwidth*iheight;m++)
    {
        if(label[m]<0)		{
			boundary[m]=1;
			R[m] = 255;
			G[m] = 0;
			B[m] = 0;
		}
        else {
			boundary[m]=0;
        }
    }

 	ProcEdgePixs(label,R,G,B,seedPos,iwidth,iheight);
 	SmootBoundaries(label,iwidth,iheight);
    
}

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{

	int iHeight;
	int iWidth;
    int i,j;

	double *R = mxGetPr(prhs[0]);
	double *G = mxGetPr(prhs[1]);
	double *B = mxGetPr(prhs[2]);
    
    double WinSize = mxGetScalar(prhs[4]);
    double lambda = mxGetScalar(prhs[3]);
 

    iWidth = mxGetM(prhs[0]);
	iHeight = mxGetN(prhs[0]);
    
	plhs[0] = mxCreateDoubleMatrix(iWidth,iHeight, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(iWidth,iHeight, mxREAL);

	double *boundary = mxGetPr(plhs[0]);
    double *label = mxGetPr(plhs[1]);

    mexPrintf("width=%d\nheight=%d\n",iWidth, iHeight);
    
    mexPrintf("WinSize = %.5f \n",WinSize);
    
    int *_label = new int[iWidth*iHeight];
    unsigned char *_boundary = new unsigned char[iWidth*iHeight];
    
	SCoW_rgb(R,G,B,_label,_boundary,iWidth,iHeight,sqrt(WinSize),lambda);
    
    for(i=0;i<iWidth*iHeight;i++)   {
        label[i] = _label[i];
        boundary[i] = _boundary[i];
    }
    
    delete []_label; _label = NULL;
    delete []_boundary; _boundary = NULL;

  }

// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}
FMatrix<float,3,3> computeTempF(vector<Match>& matches,vector<int>& inds,bool Refine=false){
    FMatrix<float,3,3> N(0.f);
    N(0,0) = 0.001;
    N(1,1) = 0.001;
    N(2,2) = 1;
    FVector<float,3> p1, p2; 
  // Construct A matrix
    Matrix<float> A;
  if (Refine){
    A=Matrix<float>(inds.size(),9);}
  else{
   A=Matrix<float>(9,9);
    //Add constraint equation
  for (int i=0; i<=8; i++){
            A(8, i) = 0;
        }}
  for(int i=0; i<inds.size(); i++) {
    p1[0] = matches[inds[i]].x1; p1[1] = matches[inds[i]].y1; p1[2]=1; 
    p2[0] = matches[inds[i]].x2; p2[1] = matches[inds[i]].y2;p2[2]=1;
    p1 = N*p1;
    p2 = N*p2;
    A(i,0) = p1[0] * p2[0];
    A(i,1) = p1[0] * p2[1];
    A(i,2) = p1[0];
    A(i,3) = p1[1] * p2[0];
    A(i,4) = p1[1] * p2[1]; 
    A(i,5) = p1[1];
    A(i,6) = p2[0];
    A(i,7) = p2[1];
    A(i,8) = 1;
  }

   
  
  Vector<float> S;       
  Matrix<float> U,Vt;
  svd(A,U,S,Vt);

  Vector<float> f = Vt.getRow(8); //The fundamental matrix is the last row of Vt

  // Construct F 
  FMatrix<float,3,3> F;
  F(0,0) = f[0]; F(0,1) = f[1]; F(0,2) = f[2];
  F(1,0) = f[3]; F(1,1) = f[4]; F(1,2) = f[5]; 
  F(2,0) = f[6]; F(2,1) = f[7]; F(2,2) = f[8];

  FVector<float,3> S2;                   
  FMatrix<float,3,3> U2, Vt2;
  svd(F,U2,S2,Vt2);

  S2[2] = 0;

  F = transpose(N)* U2*Diagonal(S2)*Vt2 *N;

  return F;

}
float computeDist(float x1,float y1,float x2,float y2,FMatrix<float,3,3> F){
    FVector<float,3> x;
    x[0] = x1;x[1] = y1;x[2] = 1;

    FVector<float,3> FTx = transpose(F) * x;
    float dist= abs(FTx[0]*x2 + FTx[1]*y2 + FTx[2]) / sqrt(FTx[0]*FTx[0] + FTx[1]*FTx[1]);
    return dist;
    
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    int nInliersMax = 0;
    for(int iter=0; iter<Niter; iter++) {
    cout<<iter<<"<- iter"<<endl;

    // Select random subset of 8 matches
    vector<int> inds;
    for(int i=0; i<8; i++){
        //inds[i] = rand() % matches.size();
        int temp = rand() % matches.size();
        inds.push_back(temp);
        }

    FMatrix<float,3,3> F =computeTempF(matches,inds);

  // Count inliers
  vector<int> Inliers;
  for(int i=0; i<matches.size(); i++) {

    float dist=computeDist(matches[i].x1,matches[i].y1,matches[i].x2,matches[i].y2,F);

    if(dist < distMax) 
      Inliers.push_back(i);
  }

  if(Inliers.size() > nInliersMax) {
    nInliersMax = Inliers.size();
    bestF = F;
    bestInliers=Inliers;
    // Update estimate of Niter
    float fracInliers = Inliers.size()/(float)matches.size();
    Niter = min(Niter, (int)(log(BETA)/log(1-pow(fracInliers,8))));
    
  }

}
    //refine F
    bestF =computeTempF(matches,bestInliers,true);


    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
        // Create a homogeneous point (x, y, 1) in the first image
        DoublePoint3 point1;
        point1[0] = x;
        point1[1] = y;
        point1[2] = 1;

        FVector<float,3> line;
        if(x>=I1.width()){
            // Display the clicked point
            drawCircle(x, y, 2, RED);
            point1[0]-= I1.width();
            // Compute the epipolar line
            line = F * point1;
            //draw the line
            float x1 = 0;
            float y1= -line[2]/line[1];;
            float x2= I1.width();
            float y2= (-line[0] *I1.width() - line[2])/line[1];
            drawLine(x1,y1, x2, y2, RED);

        }
        else{
            // Display the clicked point
            drawCircle(x, y, 2, YELLOW);
            // Compute the epipolar line
            line = transpose(F) * point1;
            //draw the line
            float x1 = I1.width();
            float y1= -line[2]/line[1];
            float x2= 2*I1.width();
            float y2= (-1)*(line[2]+line[0]*I1.width())/line[1];
            drawLine(x1,y1, x2, y2, YELLOW);
        }
        

    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}

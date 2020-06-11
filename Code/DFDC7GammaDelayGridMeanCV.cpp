/// Need to update teh time curently have delay time between +/-20, need to update the speed the reaction happens
/*

// cd /home/DavidInfortunioCygwin

//cd /home/david/BioResearch

// low-priority g++ /scratch/Infor/DF7ParallelGama.cpp  /scratch/Infor/Reaction.h -o1 -o a4.out &

 // g++ DegradAndFireParallel6withOptiSean.cpp  Reaction.h -o1 -o a4.out -std =c++11 &

// low-priority ./a4.out   &
//    fp = fopen("myfile.txt","w");

 //retVal = fwrite(buffer,sizeof(buffer),0,fp);
 //cout << "When count = 0, fwrite returned " << retVal << endl;
//#include <libstdc++>/*/
#include <iostream>
#include <algorithm>
#include <fstream>
#include <random>
#include <list>
#include "Reaction.h"
#include <algorithm>
#include <vector>
#include <thread>
#include <array>
#include <time.h>
#include<cmath>

using namespace std;

double propensity(const double&, const double&, int&);
double propensitySUM(const double&, const double&);
void printServiceLine(list<Reaction>&);
void MakeGammaFiles2(const int,const double, const double, const double, const double, const double, const double);
double my_timer();


int main() {
  const double stepsize = 0.1;
  const double startVar  = 0.1;
   const double startzz =1.9 ;
    const double addzz =1;
    const double maxzz =2;
    const double maxVar =2;
    // Exponential


    // Uniform
    uniform_real_distribution<double> distribution_uni(0.0,1.0);

    ///Commented out because not in use for Sean's experiments
    // Bernoulli

	bernoulli_distribution distribution(0.5);
   // bernoulli_distribution <double>  distribution_bern(bias);



    /// MAIN CODE SECTION

    cout << "-------------------\n-------------------\n";
        //TODO: 'Clean' variance and array size.



//https://www.gnu.org/software/parallel/
//template< class InputIt, class UnaryFunction >
//UnaryFunction for_each( InputIt first, InputIt last, UnaryFunction f );repalce the for loop below with a function and for_each
//
//const size_t nthreads = std::thread::hardware_concurrency();

  for (int j = 1; j<40; ++j) {
    MakeGammaFiles2(j,stepsize,startVar, maxVar, startzz, addzz, maxzz  );  }


    cout << "-------------------\n DONE\n-------------------\n";



    return 0;
}



/*double propensity(const double& proteinCount, const double& lambda, int i){

    switch (i){
        case 1:
            return lambda;
        case 2:
            return lambda * proteinCount;
        default:
            cout << "Not a propensity function." << endl;
            return -999;
    }
}*/

void printServiceLine(list<Reaction>& serviceLine ) {
    if (serviceLine.empty()) {
        cout << "No reactions in service line." << endl << endl;
    }
    else{
        list<Reaction>::iterator it;
        int count = 1;
        for (it = serviceLine.begin(); it != serviceLine.end(); ++it) {
            cout << "Scheduled Reaction " << count << ": " << endl
                 << "Scheduled Completion Time: " << it->getCompleationTime() << endl
                 << "Scheduled Reaction Type: " << it->getReactionType() << endl << endl;
            count++;
        }
    }


}
// [[Rcpp::export]]
void MakeGammaFiles2(const int j, const double startVar, const double stepsize,  const double maxVar,  const double startzz, const double addzz, const double maxzz  ) {

//#include <libstdc++>
#include <iostream>
#include <algorithm>

#include <fstream>
#include <random>
#include <list>
#include "Reaction.h"

  double propensity(const double&, const double&, int&);
  double propensitySUM(const double&, const double&);
  void printServiceLine(list<Reaction>&);
// double skew = 1; // This defines the stepped increase on the 'thirdMoment'. MUST divide 4 evenly. Input
 // int number_of_SSAD = 1;

 // Input This defines the number of times the SSA-D is run with the set parameters to obtain an average system time. DMI lets make this less than 10000

  double long  pTau1, pdeg2, pdeg1, pTot; // holds the odds different propensities
  std::string output1,output2;
  long double step,systemTime = 0; // Current System Time.
  int proteinCount = 0; // Current Protein Count.
 // int oldProteinCount = 0; // Current Protein Count.
  //int maxProtein; // Number of proteins that trigger a signal.
// int maxCount = 5;
//int Increase = 0;
const float endTime0 =8639.9 ;//4320 is 3 days. 432 430 changed this May 9 for time series
//const float endTime0 =399.9 ;//4320 is 3 days. 432 430
const float endTime =endTime0+0.1;
  long double tau;// tau Amount of time from Current Time to next reaction arrival.  deg1 : variable degrade rate based on protein level . deg2 : constant deg rate
  long double oneOver;

  //double   deg2;
  long double minCompletionTime = 0.0000001; // Shortest time until next scheduled reaction completion.
  // double degradePropensity=.7;
long  double sigma;
  long double theta;


 // int arraySize = (9/skew) + 1;

  int maxProduce =10; // C0 Cnot, paper is 10
  const long double rootAlpha0 = 17.32050808; // =sqrt(300); // paper was 300, need to take sqrt for letter section code squaring
  long double rootAlpha = rootAlpha0;
  //const long double yr0 = 80; //paper had it as 80
  long double yr = 80;
  //long double  yr_Over_yr0 ;
  //const long double beta0 = .1;
   long double beta = .1;
  const long double rnot = 1;
  long double Galpha, Gbeta;
  list<Reaction> serviceLine; // This 'doubly-linked-list' is used to store reactions in the service line


 //double thirdMoment = 18; // Called 'alpha' in Sean's notation; and is used to set the third moment for the delay distribution.


 // const double bias = 0.5; // This sets the weight needed for weighted-coin in the bernoulli distribution.
  default_random_engine generator (j*time(NULL)); // setting the seed for the generators.
  uniform_real_distribution<double> distribution_uni(0.0,1.0);

  ///Commented out because not in use for Sean's experiments
  // Bernoulli

    bernoulli_distribution distribution(0.5);
  // bernoulli_distribution <double>  distribution_bern(bias);

  for( double zz = startzz ; zz <=maxzz; zz=zz+addzz){//for ^1.5

    //
    /*
     //First I want to hold the dimentionless numbers constant and make sure things work as dimentionless constants before I make major changes
    //The curent setup will have the rate of change in protein level off by a feactor of yr

    //yr = yr0*zz;
   // yr_Over_yr0 = yr0/(yr);

    //rootAlpha = rootAlpha0*  sqrt(yr_Over_yr0);

    //beta = beta0*  yr_Over_yr0;
     */

    for (double i = startVar; i <= maxVar; i = i + stepsize) {//for ^2
      step = i;
      ofstream outfile;
      output1 = "/tmp/Infort/Gamma2/" ; //directory
      //output1 = "/tmp/Infort/GamaT/" ; //directory
      Galpha = 1/(zz*zz);
       Gbeta = step*(zz*zz);
	    //delay = Galpha*Gbeta
	    //cv = 1/sqrt(galpha)

       //
       // For the following setup, having a zz value greater than 3 is too hight and 1 is too low but 2 not a perfect fit
       // Galpha =16/pow(step,(zz-1));
       // Gbeta = pow(step,zz);


       //adding 1/steps to galpha decreases the time by 2 and the std by sqrt2 therby reducing the height with each step
       //adding steps to gbeta increases the time by 2 and the std by 2 therby increasing the height with each step
      //Gbeta = step*step*zz;
     // Gbeta = step*step*step*zz;
      //mean time = Gbeta*Galpha
      // std = Gbeta *sqrt(Galpha)
      //want to 4 x std and 2 x time
      // time^2 = std
      // (Gbeta*Galpha)^2 = cc*Gbeta *sqrt(Galpha)
      // Gbeta *sqrt(Galpha^3) = cc = 1/zz


      output2 = "zz"+to_string(zz)+"t"+to_string(Galpha*Gbeta)+"DelTime"+to_string(Galpha*Gbeta*Gbeta)+ "=Variance," + to_string(j)+ "_Gamma," +to_string(1/sqrt(Galpha))+"=CV.csv"; //stepsize
         //output2 = to_string(step)+ "=variance," + to_string(j)+ "_Gama," +to_string(zz)+"_yr_multiple.csv"; //stepsize

      output1 = output1 + output2 ;
      outfile.open(output1);
      proteinCount = 100;
     // proteinCount = 100*zz; // might have to change this when playing with dimentionless numbers
      systemTime = 0;
      serviceLine.clear();
      pdeg2  = pdeg1 = 0.0;


      outfile << "Time ,Proteins,peak "<<  endl;//Checking every gama distribution



      //double skew = 1; // This defines the stepped increase on the 'thirdMoment'. MUST divide 4 evenly. Input
      // int number_of_SSAD = 1;// Input This defines the number of times the SSA-D is run with the set parameters to obtain an average system time. DMI lets make this less than 10000
           gamma_distribution<double> delayGAMMA(Galpha,Gbeta);// gama(a,b) has mean a*b, variance =a*b*b, std = b*sqrt(a)
      //make mean drop as yr increases to keep dimentionless numbers the same
      //want b*sqrt(a)/b*a = 1/sqrt(a)  = 1/sqrt(1/step )= sqrt(step ) and a*b = zz
      //gamma_distribution<double> delayGAMMA(1/step,step);// gama(a,b) has mean a*b, variance =a*b*b
      // gamma_distribution<double> delayGAMMA(1/step,step*(zz));// gama(a,b) has mean a*b, variance =a*b*b, std = b*sqrt(a)
      //make mean drop as yr increases to keep dimentionless numbers the same
      //want b*sqrt(a)/b*a = 1/sqrt(a)  = 1/sqrt(1/step )= sqrt(step ) and a*b = zz

      //https://www.gnu.org/software/parallel/
      //template< class InputIt, class UnaryFunction >
      //UnaryFunction for_each( InputIt first, InputIt last, UnaryFunction f );repalce the for loop below with a function and for_each
      //
      pTot = pTau1 = 1;

      while (systemTime < endTime) { // Input DMI swiched from protein number to Systime I am fine with it going over
        long double unifDist = distribution_uni(generator);


        //double u = distribution_uni(generator);


        oneOver = rootAlpha * maxProduce/(maxProduce+proteinCount);//I need to do this or there is rounding error
        pTau1 = oneOver * oneOver;
        pdeg2 = beta*proteinCount; //dilution
        pdeg1 =  yr*proteinCount/(rnot+proteinCount);//enzymatic degrad
        pTot = (pdeg2+pTau1+pdeg1);

        exponential_distribution<double> distribution_exp(pTot);
        //instead of using the r{t-tau} in the paper, use r and wait tau to produce
        //oldProteinCount <- proteinCount;
        //				if (pTau1<=.00001){//on the off chance that the value gets lost to rounding error
        //					 pTau1 = .00001;
        //					}


        tau = distribution_exp(generator);// this differes from the papers tau, the papers tau is my theta


        if ((systemTime + tau) > (minCompletionTime) && (!serviceLine.empty())) {

          systemTime = serviceLine.front().getCompleationTime(); // update System Time.
          proteinCount++; // TODO: update with state vector.
          if(systemTime < endTime){         outfile << systemTime << "," << proteinCount << ","<<theta<< endl;}



          serviceLine.pop_front(); // Removes protein from service line.

          // sets new value for next smallest service time.
          if (serviceLine.empty()) { minCompletionTime = endTime; }
          else { minCompletionTime = serviceLine.front().getCompleationTime(); }


          continue;
        }
        else{//else Legal time



          if (unifDist*pTot<pTau1  ){  //Increase in number requested
            theta =delayGAMMA(generator);// gama(a,b) has mean a*b, variance =a*b*b
            sigma = theta + systemTime + tau;



            Reaction tempReaction = Reaction(1, sigma);// more proteins here protein1
            //Reaction tempReaction2 = Reaction(2, sigma);
            // Reaction tempReaction3 = Reaction(3, sigma);
            if (serviceLine.empty()) {//if service line empty

              systemTime = systemTime + tau;
              serviceLine.push_front(tempReaction);
              minCompletionTime = serviceLine.front().getCompleationTime();

              continue;
            }//end if service line empty
            if (sigma >= serviceLine.back().getCompleationTime()) {// if back line

              systemTime = systemTime + tau;
              serviceLine.push_back(tempReaction);
              minCompletionTime = serviceLine.front().getCompleationTime();


              continue;
            }// if back line
            else { //else no back line

              list<Reaction>::iterator it;
              for (it = serviceLine.begin(); it != serviceLine.end(); ++it) { //for stuff in line
                if (sigma < it->getCompleationTime()) {
                  systemTime = systemTime + tau;
                  serviceLine.insert(it, tempReaction);


                  break;
                }
              }//end for stuff in line

            }//end else no back line
          }//end if increased number requested
          else if  (unifDist*pTot<pTau1+pdeg2  ){  //Decrease in number requested1
            proteinCount=proteinCount-1;
            systemTime = systemTime+tau;
            sigma = systemTime ;
            outfile << systemTime << "," << proteinCount <<","<<theta<<endl;






          } //end else if Decrease in number requested1
          else if  (unifDist*pTot<pTau1+pdeg2+pdeg1 ){  //Decrease in number requested
            proteinCount=proteinCount -1;
            systemTime = systemTime+tau;
            sigma = systemTime ;
            outfile << systemTime << "," << proteinCount <<","<< theta<<endl;

          }//end else if Decrease in number requested2!







        }
        minCompletionTime = serviceLine.front().getCompleationTime();


        // TODO: Delete tempReaction???

        ///}//removed? DMI



        continue;
      }//end while*/
      /// If tau IS a 'valid' choice.
      //if (distribution_bern(generator) { theta = 3 + radiusOfDelay; } // setting
      //else { theta = 3 - radiusOfDelay; }  // setting Theta


      // TODO: rewrite to record each systemTime in an array.*/
      //systemTimeAverage = systemTimeAverage + systemTime;
      //systemTimes[i][j] = systemTime;
      // outfile  << systemTime;
      //    outfile << "," << proteinCount;

      outfile.close();
    }//end outer for-loop (&2)
  }// end ^1.5

}//End Function MakeFiles
/* Borrowing start_time struct to get timer nano-resolution.
clock_getres( CLOCK_MONOTONIC, &info.start_time );
info.clock_resolution = (int)info.start_time.tv_nsec;
/ * OK. This is permanent.
clock_gettime( CLOCK_MONOTONIC, &info.start_time );

/ * Borrowing start_time struct to get timer nano-resolution.
clock_getres( CLOCK_MONOTONIC, &info.start_time );
info.clock_resolution = (int)info.start_time.tv_nsec;

clock_gettime( CLOCK_MONOTONIC, &info.start_time );
double my_timer() {
  struct timespec t;
  clock_gettime( CLOCK_MONOTONIC, &t );
  return( (double)(t.tv_sec - info.start_time.tv_sec)
            + (t.tv_nsec - info.start_time.tv_nsec)*1.0e-9 );
}// my_timer// */

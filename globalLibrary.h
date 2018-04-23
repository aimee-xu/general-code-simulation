# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <complex.h>


const long double oneOverRoot2= 0.7071067811865475244008443621048490392848;
const long double Pi = 3.14159265358979323846264338327950288419716939937510;
const long double oneoverthree=0.3333333333333333333333333333333333333333;

// this handy little function will return the value, 0 or 1, of the bit that is at the specified location
int extractBit (const int locationOfBitFromRight, const long int theEncodedNumber)
{
  return (theEncodedNumber & ( 1L << locationOfBitFromRight )) >> locationOfBitFromRight;
}

// this function "prob" will return the total probability that a given qubit WOULD BE found in n if measured
double prob (int numQ, int q, int n, double __complex__ *amp)
{	

  long int index;
  long int stateVecSize;
  stateVecSize = 1L << numQ;
  
  double myTotal=0;
  
  for (index=0; index<stateVecSize; index++){
		if (extractBit(q,index)==n)
			myTotal += creal(amp[index])*creal(amp[index]) + cimag(amp[index])*cimag(amp[index]); 
	}
	return myTotal;
}

void swap (int numQ, int q1, int q2, double __complex__ *amp)
{
    long int NumOfSquare;
    NumOfSquare = 1L<<(numQ-larger(q1,q2)-1);
    long int SizeOfSquare;
    SizeOfSquare = (1L<<larger(q1,q2))-(1L<<smaller(q1,q2));
    
    long int i; long int j;
    
    long int temp;
    temp= 1L<<smaller(q1,q2);
    for (j=0; j<NumOfSquare; j++){
        
        for (i=temp;i<SizeOfSquare+temp;i++){
            
            if (extractBit(q1,i)!=extractBit(q2,i)){
                
                double __complex__ store=*(amp+i+SizeOfSquare);
                double __complex__ store1=*(amp+i);
                
                *(amp+i)=store;
                *(amp+i+SizeOfSquare)=store1;
            }
        }
        temp+=(1L<<(larger(q1,q2)+1));
    }
}

// this applies a sigma x gate to the qth qubit
void sigma_x (int numQ, int q, double __complex__ *amp)
{	 
  long int NumOfSquare;
  NumOfSquare = 1L << (numQ-q-1);
  long int SizeOfSquare;
  SizeOfSquare = 1L << q;
  
  
  long int i; long int j;
  	
  for (j=0; j<NumOfSquare; j++){
  	  long int temp;
      temp = SizeOfSquare*2*j;
  	for (i=temp;i<SizeOfSquare+temp;i++){
  		double __complex__ store=*(amp+i+SizeOfSquare);
  		double __complex__ store1=*(amp+i);
  		
  		*(amp+i)=store;
  		*(amp+i+SizeOfSquare)=store1;
  	} 	
  }
}

// this applies a sigma y gate to the qth qubit
void sigma_y (int numQ, int q, double __complex__ *amp)
{	 
  long int NumOfSquare;
  NumOfSquare = 1L << (numQ-q-1);
  long int SizeOfSquare;
  SizeOfSquare = 1L << q;
  
  
  long int i; long int j;
  	
  for (j=0; j<NumOfSquare; j++){
  	  long int temp;
      temp = SizeOfSquare*2*j;
  	for (i=temp;i<SizeOfSquare+temp;i++){
  		double __complex__ store=-*(amp+i+SizeOfSquare)*I;
  		double __complex__ store1=*(amp+i)*I;
  		
  		*(amp+i)=store;
  		*(amp+i+SizeOfSquare)=store1;
  	} 	
  }
}

// this applies a sigma z gate to the qth qubit
void sigma_z(int numQ, int q, double __complex__ *amp)
{	 
long int stateVecSize;
stateVecSize = 1L << numQ;
  
  long int index;
  	
  for (index=0; index<stateVecSize; index++){
		 if (extractBit(q,index)==1)
		 	*(amp+index)*=(-1);
   }
}

void noisySigma_z (int numQ, int q, double gateError, double __complex__ *amp){
	
    sigma_z (numQ, q, amp);
	double ran1; double ran2;
	ran1 =((double)rand()/(double)RAND_MAX);
	ran2 =((double)rand()/(double)RAND_MAX);
    
	if (ran1<gateError){	
			if (ran2<=oneoverthree)
				sigma_x (numQ, q, amp);
			else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
				sigma_y (numQ, q, amp);
			else if (ran2>2*oneoverthree)
				sigma_z (numQ, q, amp);
	}
}

void noisySigma_y (int numQ, int q, double gateError, double __complex__ *amp){
    
    sigma_y (numQ, q, amp);
	double ran1; double ran2;
	ran1 =((double)rand()/(double)RAND_MAX);
	ran2 =((double)rand()/(double)RAND_MAX);
		
	if (ran1<gateError){	
			if (ran2<=oneoverthree)
				sigma_x (numQ, q, amp);
			else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
				sigma_y (numQ, q, amp);
			else if (ran2>2*oneoverthree)
				sigma_z (numQ, q, amp);
    }
}

void noisySigma_x (int numQ, int q, double gateError, double __complex__ *amp){
    
    sigma_x (numQ, q, amp);
	double ran1; double ran2;
	ran1 =((double)rand()/(double)RAND_MAX);
	ran2 =((double)rand()/(double)RAND_MAX);
		
	if (ran1<gateError){	
			if (ran2<=oneoverthree)
				sigma_x (numQ, q, amp);
			else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
				sigma_y (numQ, q, amp);
			else if (ran2>2*oneoverthree)
				sigma_z (numQ, q, amp);
	}
	
}

void initializationNoise (int numQ, int q, double gateError, double __complex__ *amp){
    
    double ran1; double ran2;
    ran1 =((double)rand()/(double)RAND_MAX);
    ran2 =((double)rand()/(double)RAND_MAX);
    
    if (ran1<gateError){
        if (ran2<=oneoverthree)
            sigma_x (numQ, q, amp);
        else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
            sigma_y (numQ, q, amp);
        else if (ran2>2*oneoverthree)
            sigma_z (numQ, q, amp);
    }
    
}

// this applies a hadamard gate to the qth qubit
void hadamard (int numQ, int q, double __complex__ *amp)
{	 
  long int NumOfSquare;
  NumOfSquare = 1L << (numQ-q-1);
  long int SizeOfSquare;
  SizeOfSquare = 1L << q;
  
  
  long int i; long int j;
  	
  for (j=0; j<NumOfSquare; j++){
  	  long int temp;
      temp = SizeOfSquare*2*j;
  	for (i=temp;i<SizeOfSquare+temp;i++){
  		double __complex__ store=(*(amp+i)+*(amp+i+SizeOfSquare))*oneOverRoot2;
  		double __complex__ store1=(*(amp+i)-*(amp+i+SizeOfSquare))*oneOverRoot2;
  		
  		*(amp+i)=store;
  		*(amp+i+SizeOfSquare)=store1;
  	} 	
  }
}

void noisyHadamard (int numQ, int q, double gateError, double __complex__ *amp){
    
    hadamard (numQ, q, amp);
	double ran1; double ran2;
	ran1 =((double)rand()/(double)RAND_MAX);
	ran2 =((double)rand()/(double)RAND_MAX);
	
	if (ran1<gateError){	
			if (ran2<=oneoverthree)
				sigma_x (numQ, q, amp);
			else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
				sigma_y (numQ, q, amp);
			else if (ran2>2*oneoverthree)
				sigma_z (numQ, q, amp);
	}
		
}

// this applies a cphase gate to q1 and q2
void cphase (int numQ, int q1, int q2, double __complex__ *amp)
{	
  long int stateVecSize;
  stateVecSize = 1L << numQ;
  long int index;

  for (index=0; index<stateVecSize; index++){
  	 if ((extractBit(q1,index))*(extractBit(q2,index))==1)
			*(amp+index) *= (-1);
	}
}

// this applies a cnot gate, where q1 is the control qubit, q2 is the target qubit
void cnot (int numQ, int q1, int q2, double __complex__ *amp)
{	 
  long int NumOfSquare;
  NumOfSquare = 1L << (numQ-q2-1);
  long int SizeOfSquare;
  SizeOfSquare = 1L << q2;
    
    long int i; long int j;
  	
  for (j=0; j<NumOfSquare; j++){
  	  long int temp;
      temp = SizeOfSquare*2*j;
  	for (i=temp;i<SizeOfSquare+temp;i++){
  	if (extractBit(q1,i)==1){
  		double __complex__ store=*(amp+i+SizeOfSquare);
  		double __complex__ store1=*(amp+i);
  		
  		*(amp+i)=store;
  		*(amp+i+SizeOfSquare)=store1;
  	 } 	
   }
 }
}

void noisyCnot (int numQ, int q1, int q2, double gateError, double __complex__ *amp){
	
  cnot (numQ, q1, q2, amp);
  double ran;
  ran =((double)rand()/(double)RAND_MAX);
    if (ran < gateError)
    {
        int myRnd=1 + rand()%15;
        int firstError=myRnd/4;
        int secndError=myRnd%4;
        if (firstError==1) sigma_x (numQ, q1, amp);
        if (firstError==2) sigma_y (numQ, q1, amp);
        if (firstError==3) sigma_z (numQ, q1, amp);
        if (secndError==1) sigma_x (numQ, q2, amp);
        if (secndError==2) sigma_y (numQ, q2, amp);
        if (secndError==3) sigma_z (numQ, q2, amp);
    }

}


void noisyCphase (int numQ, int q1, int q2, double gateError, double __complex__ *amp){
	
    cphase (numQ, q1,  q2, amp);
	double ran;
    ran =((double)rand()/(double)RAND_MAX);
    if (ran < gateError)
    {
        int myRnd=1 + rand()%15;
        int firstError=myRnd/4;
        int secndError=myRnd%4;
        if (firstError==1) sigma_x (numQ, q1, amp);
        if (firstError==2) sigma_y (numQ, q1, amp);
        if (firstError==3) sigma_z (numQ, q1, amp);
        if (secndError==1) sigma_x (numQ, q2, amp);
        if (secndError==2) sigma_y (numQ, q2, amp);
        if (secndError==3) sigma_z (numQ, q2, amp);
    }

}


void rotation_angle (int numQ, int q, double eta, double theta1, double theta2, double theta3, double __complex__ *amp)
{	 
  long int NumOfSquare; 
  NumOfSquare = 1L << (numQ-q-1);
  long int SizeOfSquare;
  SizeOfSquare = 1L << q;
  
  double __complex__ a; double __complex__ b; 
  double __complex__ c; double __complex__ d;
    
  theta1=(theta1+eta*(rand()/(RAND_MAX-1.0)*2-1))*Pi/180;
  theta2=(theta2+eta*(rand()/(RAND_MAX-1.0)*2-1))*Pi/180;
  theta3=(theta3+eta*(rand()/(RAND_MAX-1.0)*2-1))*Pi/180;
  
  a=cos(theta1)*(cos(theta2)+I*sin(theta2));
  d=cos(theta1)*(cos(theta2)-I*sin(theta2));
  b=sin(theta1)*(cos(theta3)+I*sin(theta3));
  c=sin(theta1)*(-cos(theta3)+I*sin(theta3));
  
  long int i; long int j;
  	
  for (j=0; j<NumOfSquare; j++){
  	  long int temp;
      temp = SizeOfSquare*2*j;
  	for (i=temp;i<SizeOfSquare+temp;i++){
  		double __complex__ store=*(amp+i)*a+*(amp+i+SizeOfSquare)*b;
  		double __complex__ store1=*(amp+i)*c+*(amp+i+SizeOfSquare)*d;
  		
  		*(amp+i)=store;
  		*(amp+i+SizeOfSquare)=store1;
  	} 	
  }
}

//this function returns the collapsed state if we want the qth qubit to be n(0 or 1) after the measurement
void measure(int numQ, int q, int n, double __complex__ *amp)
{
    long int stateVecSize;
  // dimension of the state vector
  stateVecSize = 1L << numQ;
   
  long int index;
  long double pro= prob(numQ, n, q, amp);
  if (pro<=1e-20){ printf("very very small probability:%Lf", pro); exit(1); }
  long double recipnorm=1/sqrt(pro);
  	
  for (index=0; index<stateVecSize; index++){
  	if (extractBit(q,index)==n){
        *(amp+index)=(*(amp+index))*recipnorm;  
  	}else{
  		*(amp+index)=0;
  	} 		
  }	
}

//this function returns the collapsed state if we want the qth qubit to be absolutely n(0 or 1) after the measurement
void measureAbsolute(int numQ, int q, int n, double __complex__ *amp)
{
    long int stateVecSize;
  // dimension of the state vector
  stateVecSize = 1L << numQ;
   
  long int index;
  long double pro= prob(numQ, q, n, amp);
  if (pro<=1e-20){ sigma_x(numQ, q, amp);}
  else
  {	  long double recipnorm=1/sqrt(pro);
		
	  for (index=0; index<stateVecSize; index++){
		if (extractBit(q,index)==n){
			*(amp+index)=(*(amp+index))*recipnorm;  
		}else{
			*(amp+index)=0;
		} 		
	  }	
  }
}

void noisyMeasure(int numQ, int q, int n, double measureError, double __complex__ *amp){
	
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    if (ran<measureError)
        sigma_x (numQ, q, amp);
				
	measure(numQ, q, n,amp);
}

int randomCollapse(int numQ, int q, double measureError, double __complex__ *amp)
{
    
  int selectedOutcome = -1;
  long double pro= prob(numQ, q, 0, amp);
  if (pro<=1e-15) selectedOutcome=1; 		//these two special cases should not really be necessary but they catch situations that
  if ((1-pro)<=1e-15) selectedOutcome=0;		//are likely to be associated with small numerical errors where the 'correct' prob is 0
   if (selectedOutcome==-1){
   		if (pro>((double)rand()/(double)RAND_MAX)){
   			selectedOutcome=0;
   		}else{
   			selectedOutcome=1;
   		}
   }
    
    measureAbsolute(numQ, q, selectedOutcome ,amp);
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    if (ran<measureError){
        if (selectedOutcome==0)
            selectedOutcome=1;
        else
            selectedOutcome=0;
        }

  return selectedOutcome;
}

void dephasingModel (int numQ, int q, double t,  double __complex__ *amp)
{
    double pd = 0.5*(1-exp(-t));
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    
    if (ran<pd)
        sigma_z (numQ, q, amp);
}

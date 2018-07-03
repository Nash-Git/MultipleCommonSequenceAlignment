#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include<ctime>
#include <algorithm>
using namespace std;

#define MATCH      4
#define MISMATCH  -1
#define INDEL     -3

vector<char> align_seq1;
vector<char> align_seq2;

// This fuction calculate match and mismatch for global alignment

int match_mismatch(char a,char b){

  int result;
  if(a==b){
      result=MATCH;
    }
  else{
      result=MISMATCH;
    }
  return result;
}
// This fucntion return max value for recurence relation

int max_S(int a[],int length){

  int max = a[0];       

  for(int i = 1; i<length; i++)
  {
      if(a[i] > max)
	   {
		max = a[i];
	   }
  }
  return max;                    
}

// Retuen min value for selecting centre star sequence
int min_S(int a[],int length){

  int min = 100000;
  int l;

  for(int i = 0; i<length; i++)
  {
      if(a[i] < min)
	   {
		min = a[i];
        l=i;
	   }
  }
 
  return l;                    
}
// For printing every alignment with centre star sequence
void print_alignment (vector <char> align_seq1, vector <char> align_seq2)
{

	for(int i=align_seq1.size()-1;i>=0;i--)
		cout<<align_seq1[i];
		cout<<endl;
	for(int i=align_seq2.size()-1;i>=0;i--)
		cout<<align_seq2[i];
		cout<<endl;
		cout<<"Happy"<<endl;
}

// backtracking the global alignment
void backtrack( vector<char> Seq1, vector<char> Seq2, vector< vector <int> > S,int i, int j)
{
	int k=0;
	int score,count=0;
	
	align_seq1.clear();
	align_seq2.clear();

	while (i > 0 || j > 0)
   {
		if (i > 0 && j > 0) 
		{
			score = S[i-1][j-1];
			if (Seq1[i-1] == Seq2[j-1])
				score += MATCH; 
			else 
				score += MISMATCH;
			if (score == S[i][j]) 
			{
				align_seq1.insert(align_seq1.begin()+k,Seq1[i-1]);
				align_seq2.insert(align_seq2.begin()+k,Seq2[j-1]);
				k++; 
				i--; 
				j--;
				continue;
			 }
		   }

		if(i>0)
		{
			if (S[i-1][j] + INDEL == S[i][j]) 
			{
			
				align_seq1.insert(align_seq1.begin()+k,Seq1[i-1]);		
				align_seq2.insert(align_seq2.begin()+k,'-');
				k++; 
				i--; 
				continue;
			}
		 }

		if(j>0)
		{
			if (S[i][j-1] + INDEL == S[i][j]) 
			{
				align_seq1.insert(align_seq1.begin()+k,'-');
				align_seq2.insert(align_seq2.begin()+k,Seq2[j-1]);
				k++;
				j--;
				continue;
			}
		}
	}
	
	print_alignment(align_seq1,align_seq2);
	
}

// calculating matrix for global alignment
void globalmatrix(const vector<char> &Seq1,const vector<char> Seq2)
{   
        vector< vector <int> >  S (Seq1.size()+1 , vector<int> (Seq2.size() +1 ) );
		int temp[3];
		int sc;
		int len_seq1 = Seq1.size();
		int len_seq2 = Seq2.size();
		
		S[0][0] = 0;
		for (int i = 1; i <= len_seq1; i++)  
			S[i][0] = i * INDEL;
		for(int j = 1; j <= len_seq2; j++)  
			S[0][j] = j * INDEL;

        for(int i =1; i<=len_seq1; i++)
		{
			for( int j =1; j<=len_seq2; j++)
			{
				temp[0]=  S[i-1][j-1]+match_mismatch(Seq1[i-1],Seq2[j-1]); // add match and mismatch
				temp[1] = S[i-1][j] + INDEL;               
				temp[2] = S[i][j-1] + INDEL; 
				S[i][j] = max_S(temp,3);
								
			 }
		
		}
		backtrack(Seq1,Seq2,S,len_seq1,len_seq2);
  
}
// Print the final alignment after merging
void print(vector< vector <char> > final){

		for (int i=0;i<24;i++)
		{
			for (int j=0;j<final[i].size();j++)
			{
			cout<<final[i][j];
			}
			cout<<endl;
		}

}

// Selecting star,adding indel, and merging
void star (int c[][30],int n,vector< vector <char> > input){
		int d[30];
		int i,j;
		vector < vector<char> >  final(align_seq1.size()+35 , vector<char> (align_seq2.size() +1 ));
		int k=0;
		vector <char> check_indel;
		for ( i=0;i<n;i++)
		{   
			int sum=0;
			for ( j=0;j<n;j++)
			{
				sum=sum+c[i][j]; // Calculate the sume to get the lowest mismatch
						
			}
			d[i]=sum;
			
		}

		int min_d=min_S(d,n); // Select the star sequence
			
		check_indel=input[min_d];  
		int track,x;
		for ( i=0;i<n;i++)
		{  
			
			if(i!=min_d){
				if (input[min_d].size()>=align_seq1.size()) 	{
							
				globalmatrix(input[min_d],input[i]);
				reverse(align_seq2.begin(),align_seq2.end());
				final[k]=align_seq2;
				k++;
				}
				else 
				{
					
					
					
					if (check_indel.size()<align_seq1.size())
					{
							
							for (int d=k-1;d>=0;d--)
							{
								reverse(align_seq1.begin(),align_seq1.end());
					  			globalmatrix(align_seq1,input[d]);
								reverse(align_seq2.begin(),align_seq2.end());
								final[d]=align_seq2;
							}
						
					}
					check_indel=align_seq1;
					reverse(align_seq1.begin(),align_seq1.end());
					globalmatrix(align_seq1,input[i]);
					reverse(align_seq2.begin(),align_seq2.end());
					final[k]=align_seq2;
					k++;

					
				}
			}
			
		
		}
		
		reverse(align_seq1.begin(),align_seq1.end());
		final[k]=align_seq1;
		final[k-1]=align_seq2;
		if (check_indel.size()<align_seq1.size()){
						for (int d=k-2;d>=0;d--)
							{
								
					  			globalmatrix(align_seq1,input[d]);
								reverse(align_seq2.begin(),align_seq2.end());
								final[d]=align_seq2;
								reverse(align_seq1.begin(),align_seq1.end());
							}
		
		}
				
print(final);
}
  
// Initial mismatch calculation for matrix creation
int centrematrix(const vector<char> &Seq1,const vector<char> Seq2)
{   
        int count=0,x,y,s;
		x=Seq1.size();
		y=Seq2.size();
		if (x<=y)
			s=x;
		else
			s=y;

		for (int l=0;l<s;l++)
		{
	
			if (Seq1[l]!=Seq2[l])
			count++;
		}
		return count;
}


// Read 24 input file sequence
void readfile()
{
	 int i,j,n;
	 n=24;
	 int c[30][30];
	 
	 ifstream one("C://test/f");
     istream_iterator<char> start(one), end;
     vector<char> x1(start, end);
	 ifstream two("C://test/f1");
     istream_iterator<char> start1(two), end1;
     vector<char> x2(start1, end1);

	 ifstream three("C://test/f2");
     istream_iterator<char> start2(three), end2;
     vector<char> x3(start2, end2);
	 ifstream four("C://test/f3");
     istream_iterator<char> start3(four), end3;
     vector<char> x4(start3, end3);
	 ifstream five("C://test/f4");
     istream_iterator<char> start4(five), end4;
     vector<char> x5(start4, end4);

	 ifstream six("C://test/f5");
     istream_iterator<char> start5(six), end5;
     vector<char> x6(start5, end5);
	 ifstream seven("C://test/f6");
     istream_iterator<char> start6(seven), end6;
     vector<char> x7(start6, end6);
	 ifstream eight("C://test/f7");
     istream_iterator<char> start7(eight), end7;
     vector<char> x8(start7, end7);
	 
	 ifstream nine("C://test/f8");
     istream_iterator<char> start8(nine), end8;
     vector<char> x9(start8, end8);
	 
	 ifstream ten("C://test/f9");
     istream_iterator<char> start9(ten), end9;
     vector<char> x10(start9, end9);

	 ifstream eleven("C://test/f10");
     istream_iterator<char> start10(eleven), end10;
     vector<char> x11(start10, end10);

	 ifstream twelve("C://test/f11");
     istream_iterator<char> start11(twelve), end11;
     vector<char> x12(start11, end11);

	 ifstream thirteen("C://test/f12");
     istream_iterator<char> start12(thirteen), end12;
     vector<char> x13(start12, end12);

	 ifstream fourteen("C://test/f13");
     istream_iterator<char> start13(fourteen), end13;
     vector<char> x14(start13, end13);

	 ifstream fifteen("C://test/f14");
     istream_iterator<char> start14(fifteen), end14;
     vector<char> x15(start14, end14);

	 ifstream sixteen("C://test/f15");
     istream_iterator<char> start15(sixteen), end15;
     vector<char> x16(start15, end15);

	 ifstream seventeen("C://test/f16");
     istream_iterator<char> start16(seventeen), end16;
     vector<char> x17(start16, end16);
	 
	 ifstream eighteen("C://test/f17");
     istream_iterator<char> start17(eighteen), end17;
     vector<char> x18(start17, end17);

	 ifstream nineteen("C://test/f18");
     istream_iterator<char> start18(nineteen), end18;
     vector<char> x19(start18, end18);

	 ifstream twenty("C://test/f19");
     istream_iterator<char> start19(twenty), end19;
     vector<char> x20(start19, end19);

	 ifstream twenty_one("C://test/f20");
     istream_iterator<char> start20(twenty_one), end20;
     vector<char> x21(start20, end20);

	 ifstream twenty_two("C://test/f21");
     istream_iterator<char> start21(twenty_two), end21;
     vector<char> x22(start21, end21);

	 ifstream twenty_three("C://test/f22");
     istream_iterator<char> start22(twenty_three), end22;
     vector<char> x23(start22, end22);

	 ifstream twenty_four("C://test/f23");
     istream_iterator<char> start23(twenty_four), end23;
     vector<char> x24(start23, end23);

	 vector < vector<char> >  input;

	 // Pushing all the input of file to a 2D vector input    

	for(int i=0;i<24;i++)
	{
	   input.push_back(vector<char> ());
	   if (i==0)
	  {
	   for(j=0;j<x1.size();j++)
	   {
		   input[i].push_back(x1[j]);
	       
	   }
	  }

	   else if (i==1)
	   {
	   
	   for( j=0;j<x2.size();j++)
	   {
		   input[i].push_back(x2[j]);
	       
	   }
	   
	   }

	   else if (i==2)
	   {
	   
	   for( j=0;j<x3.size();j++)
	   {
		   input[i].push_back(x3[j]);
	       
	   }
	   
	   }
	   else if (i==3)
	   {
	   
	   for( j=0;j<x4.size();j++)
	   {
		   input[i].push_back(x4[j]);
	       
	   }
	   
	   }
	   else if (i==4)
	   {
	   
	   for( j=0;j<x5.size();j++)
	   {
		   input[i].push_back(x5[j]);
	       
	   }
	   
	   }
	   else if (i==5)
	   {
	   
	   for( j=0;j<x6.size();j++)
	   {
		   input[i].push_back(x6[j]);
	       
	   }
	   
	   }
	   else if (i==6)
	   {
	   
	   for( j=0;j<x7.size();j++)
	   {
		   input[i].push_back(x7[j]);
	       
	   }
	   
	   }
	   else if (i==7)
	   {
	   
	   for( j=0;j<x8.size();j++)
	   {
		   input[i].push_back(x8[j]);
	       
	   }
	   
	   }

	   else if (i==8)
	   {
	   
	   for( j=0;j<x9.size();j++)
	   {
		   input[i].push_back(x9[j]);
	       
	   }
	   
	   }

	   else if (i==9)
	   {
	   
	   for( j=0;j<x10.size();j++)
	   {
		   input[i].push_back(x10[j]);
	       
	   }
	   
	   }

	   else if (i==10)
	   {
	   
	   for( j=0;j<x11.size();j++)
	   {
		   input[i].push_back(x11[j]);
	       
	   }
	   
	   }

	   else if (i==11)
	   {
	   
	   for( j=0;j<x12.size();j++)
	   {
		   input[i].push_back(x12[j]);
	       
	   }
	   
	   }

	   else if (i==12)
	   {
	   
	   for( j=0;j<x13.size();j++)
	   {
		   input[i].push_back(x13[j]);
	       
	   }
	   
	   }

	   else if (i==13)
	   {
	   
	   for( j=0;j<x14.size();j++)
	   {
		   input[i].push_back(x14[j]);
	       
	   }
	   
	   }

	   else if (i==14)
	   {
	   
	   for( j=0;j<x15.size();j++)
	   {
		   input[i].push_back(x15[j]);
	       
	   }
	   
	   }

	   else if (i==15)
	   {
	   
	   for( j=0;j<x16.size();j++)
	   {
		   input[i].push_back(x16[j]);
	       
	   }
	   
	   }

	   else if (i==16)
	   {
	   
	   for( j=0;j<x17.size();j++)
	   {
		   input[i].push_back(x17[j]);
	       
	   }
	   
	   }

		else if (i==17)
	   {
	   
	   for( j=0;j<x18.size();j++)
	   {
		   input[i].push_back(x18[j]);
	       
	   }
	   
	   }

	   else if (i==18)
	   {
	   
	   for( j=0;j<x19.size();j++)
	   {
		   input[i].push_back(x19[j]);
	       
	   }
	   
	   }

		else if (i==19)
	   {
	   
	   for( j=0;j<x20.size();j++)
	   {
		   input[i].push_back(x20[j]);
	       
	   }
	   
	   }

	   else if (i==20)
	   {
	   
	   for( j=0;j<x21.size();j++)
	   {
		   input[i].push_back(x21[j]);
	       
	   }
	   
	   }

		else if (i==21)
	   {
	   
	   for( j=0;j<x22.size();j++)
	   {
		   input[i].push_back(x22[j]);
	       
	   }
	   
	   }

	   else if (i==22)
	   {
	   
	   for( j=0;j<x23.size();j++)
	   {
		   input[i].push_back(x23[j]);
	       
	   }
	   
	   }

		else if (i==23)
	   {
	   
	   for( j=0;j<x24.size();j++)
	   {
		   input[i].push_back(x24[j]);
	       
	   }
	   
	   }

	}
 // generate the centre stat matrix
	for ( i=0;i<n;i++)
	{
		for ( j=i;j<n-1;j++)
		{
			if (i==j){
				c[i][j]=0;
				c[i][j+1]=centrematrix(input[i],input[j+1]);
				c[j+1][i]=c[i][j+1];
					 }
			else 
				c[i][j+1]=centrematrix(input[i],input[j+1]);
				c[j+1][i]=c[i][j+1];
				 
		}
		
		
		if (i==(n-1))
		c[i][j]=0;
	}

	// Printing the matrix
	  for ( i=0;i<n;i++)
	{
		for ( j=0;j<n;j++)
		{
			cout<<c[i][j]<<"\t";
		}
		cout<<endl;
	}

 star(c,n,input); // calling the stat method with matrix, number of input and all seequence which stored in input
	
}

int main()
{
	 clock_t startTime = clock();

	 readfile();
	  
	 
 	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;

	 system("pause");
	return 0;
}
  
library(Rcpp)
cppFunction("
// C++ program for quantile normalization
std::vector<std::vector<double> >  rowmeanquantile(std::vector<std::vector<double> > mat);
std::vector<std::vector<double> > transpose(	std::vector<std::vector<double> > mat);
std::vector<std::vector<double> >  sortCol(std::vector<std::vector<double> > mat);
void RowWiseSort(std::vector<std::vector<double> >& B);


std::vector<std::vector<double> >  rowmeanquantile(std::vector<std::vector<double> > mat)
            {double sum;
             int N = mat.size();//row
	int M = mat[0].size();//column
            double mean_val[N];     
              
 std:: map<double, double> mean;  
  std::map<double,double>::iterator it;
 std::vector<std::vector<double> > df= mat;
  df=sortCol(mat);
  
            for (int i = 0; i < N; i++) {
                    sum=0;
                    for (int j = 0; j < M; j++) {
                        sum+=df[i][j];
                    }
                   mean_val[i]=sum/M;
                      
                    }
               
            
            for (int i = 0; i < N; i++) {
                    for (int j = 0; j < M; j++) {
                         mean.insert(std::pair<double, double>(df[i][j],mean_val[i]));     
                         
                    }}
        for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
            for ( auto it = mean.begin(); it != mean.end(); ++it  )
            {
         if(mat[i][j]==it->first)
         {
             mat[i][j]=it->second;
         }
            }
        
		
        }
       
		}

            return mat;
            }

// Function to find the transpose
// of the matrix mat[]
std::vector<std::vector<double> > transpose(
	std::vector<std::vector<double> > mat)
{
     int row = mat.size();//row
	int col= mat[0].size();//column

	// Stores the transpose
	// of matrix mat[][]
	std::vector<std::vector<double> > tr(
		col, std::vector<double>(row));

	// Traverse each row of the matrix
	for (int i = 0; i < row; i++) {

		// Traverse each column of the matrix
		for (int j = 0; j < col; j++) {

			// Transpose matrix elements
			tr[j][i] = mat[i][j];
		}
	}

	return tr;
}

// Function to sort the given
// matrix in row wise manner
void RowWiseSort(std::vector<std::vector<double> >& B)
{
	// Traverse the row
	for (int i = 0; i < (int)B.size(); i++) {

		// Row - Wise Sorting
		sort(B[i].begin(), B[i].end());
	}
}

// Function to print the matrix
// in column wise sorted manner
std::vector<std::vector<double> >  sortCol(std::vector<std::vector<double> > mat)
{

	// Function call to find transpose
	// of the the matrix mat[][]
	std::vector<std::vector<double> > B
		= transpose(mat);

	// Sorting the matrix row-wise
	RowWiseSort(B);

	// Calculate transpose of B[][]
	mat = transpose(B);

	return mat;
}

")



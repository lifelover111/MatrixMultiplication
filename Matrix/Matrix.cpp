#include <omp.h>
#include <stdio.h>
#include <conio.h>
#include <chrono>
#include <random>
#include <iostream>


using namespace std;

double** MultiplyMatrix(double** A, double** B, int size);
void DisplayMatrix(double** A, int size);
double** GenerateMatrix(int size);
double*** SplitMatrix(double**& matrix, int& size);
double** CollectMatrix(double*** blocks, int size);
double** Sum(double** A, double** B, int size);
double** Subtract(double** A, double** B, int size);
double** Strassen(double** A, double** B, int& size);
double** ParallelMultiply(double** A, double** B, int size);
double** ParallelStrassen(double** A, double** B, int& size);
double*** ParallelSplit(double**& matrix, int& size);
double** ParallelCollect(double*** blocks, int size);

int main()
{
    int size, mode;
    double** A,** B,** C;
    double t1, t2;
    while (true)
    {

        cout << "Enter size of matrix(size x size): ";
        cin >> size;
        A = GenerateMatrix(size);
        B = GenerateMatrix(size);
        cout << "Random matrix are generated \n";
        cout << "Enter mode (1 - common multiply; 2 - parallel multiply; 3 - Strassen; 4 - parallel Strassen): ";
        cin >> mode;
        switch (mode)
        {
        case 1:
            t1 = omp_get_wtime();
            C = MultiplyMatrix(A, B, size);
            t2 = omp_get_wtime();
            cout << "The multiply time: " << t2 - t1 << "\n\n";
            break;
        case 2:
            t1 = omp_get_wtime();
            C = ParallelMultiply(A, B, size);
            t2 = omp_get_wtime();
            cout << "The multiply time: " << t2 - t1 << "\n\n";
            break;
        case 3:
            t1 = omp_get_wtime();
            C = Strassen(A, B, size);
            t2 = omp_get_wtime();
            cout << "The multiply time: " << t2 - t1 << "\n\n";
            break;
        case 4:
            t1 = omp_get_wtime();
            C = ParallelStrassen(A, B, size);
            t2 = omp_get_wtime();
            cout << "The multiply time: " << t2 - t1 << "\n\n";
            break;
        default:
            break;
        }

        delete[] A, B, C;

    }
}

double** MultiplyMatrix(double** A, double** B, int size)
{
    double** result = new double* [size];
    for (int i = 0; i < size; i++)
    {
        result[i] = new double[size];
        for (int j = 0; j < size; j++)
        {
            result[i][j] = 0;
        }
    }

    for (int r = 0; r < size; r++)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                result[i][j] += A[i][r] * B[r][j];
            }
        }
    }

    return result;
}

void DisplayMatrix(double** A, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

double** GenerateMatrix(int size)
{
    random_device rd;
    mt19937::result_type seed = rd() ^ ((mt19937::result_type)
        chrono::duration_cast<chrono::seconds>(chrono::system_clock::now().time_since_epoch()).count() + (mt19937::result_type)
        chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
    mt19937 gen(seed);
    uniform_real_distribution<> distrib(0, size);


    double** M = new double* [size];
    for (int i = 0; i < size; i++)
    {
        M[i] = new double[size];
        for (int j = 0; j < size; j++)
        {
            M[i][j] = distrib(gen);
        }
    }
    return M;
}

double*** SplitMatrix(double**& matrix, int& size)
{
    if (size % 2 == 1)
    {
        double** temp = new double* [size + 1];
        for (int i = 0; i < size + 1; i++)
        {
            temp[i] = new double[size + 1];
            for (int j = 0; j < size + 1; j++)
            {
                if ((i == size) || (j == size))
                {
                    temp[i][j] = 0;
                }
                else
                {
                    temp[i][j] = matrix[i][j];
                }
            }
        }
        matrix = temp;
        size = size + 1;
    }


    double*** blocks = new double** [4];

    for (int k = 0; k < 4; k++)
    {
        blocks[k] = new double* [size/2];
        for (int i = 0; i < size/2; i++)
        {
            blocks[k][i] = new double[size/2];
            for (int j = 0; j < size/2; j++)
            {
                blocks[k][i][j] = matrix[(size / 2) * (k / 2) + i][(size / 2) * (k % 2) + j];
            }
        }
    }

    return blocks;
}

double** CollectMatrix(double*** blocks, int size)
{
    double** matrix = new double* [size];
    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
    }

    for (int k = 0; k < 4; k++)
    {
        for (int i = 0; i < size / 2; i++)
        {
            for (int j = 0; j < size / 2; j++)
            {
                matrix[(size / 2) * (k / 2) + i][(size / 2) * (k % 2) + j] = blocks[k][i][j];
            }
        }
    }

    return matrix;
}


double** Sum(double** A, double** B, int size)
{
    double** result = new double* [size];
    for (int i = 0; i < size; i++)
    {
        result[i] = new double[size];
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}

double** Subtract(double** A, double** B, int size)
{
    double** result = new double* [size];
    for (int i = 0; i < size; i++)
    {
        result[i] = new double[size];
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}


double** Strassen(double** A, double** B, int& size)
{
    int size1 = size;
    double*** blocksA = SplitMatrix(A, size);
    double*** blocksB = SplitMatrix(B, size1);
    
    double** D = MultiplyMatrix(Sum(blocksA[0], blocksA[3], size / 2), Sum(blocksB[0], blocksB[3], size / 2), size / 2);
    double** H2 = MultiplyMatrix(Sum(blocksA[2], blocksA[3], size / 2), blocksB[0], size / 2);
    double** V2 = MultiplyMatrix(blocksA[0],Subtract(blocksB[1], blocksB[3], size/2), size / 2);
    double** V1 = MultiplyMatrix(blocksA[3], Subtract(blocksB[2], blocksB[0], size / 2), size / 2);
    double** H1 = MultiplyMatrix(Sum(blocksA[0], blocksA[1], size / 2), blocksB[3], size / 2);
    double** D2 = MultiplyMatrix(Subtract(blocksA[2], blocksA[0], size / 2), Sum(blocksB[0], blocksB[1], size / 2), size / 2);
    double** D1 = MultiplyMatrix(Subtract(blocksA[1], blocksA[3], size / 2), Sum(blocksB[2], blocksB[3], size / 2), size / 2);


    double*** R = new double** [4];
    R[0] = Subtract(Sum(Sum(D, D1, size / 2), V1, size / 2), H1, size / 2);
    R[1] = Sum(V2, H1, size / 2);
    R[2] = Sum(V1, H2, size / 2);
    R[3] = Subtract(Sum(Sum(D, D2, size / 2), V2, size / 2), H2, size / 2);
    double** result = CollectMatrix(R, size);
    delete[] D, D1, D2, V1, V2, H1, H2, R[0], R[1], R[2], R[3];
    return result;
}


double** ParallelMultiply(double** A, double** B, int size)
{
    double** result = new double* [size];
    for (int i = 0; i < size; i++)
    {
        result[i] = new double[size];
        for (int j = 0; j < size; j++)
        {
            result[i][j] = 0;
        }
    }
#pragma omp parallel for
        for (int r = 0; r < size; r++)
        {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    result[i][j] += A[i][r] * B[r][j];
                }
            }
        }
    


    return result;
}

double** ParallelStrassen(double** A, double** B, int& size)
{
    int size1 = size;
    double*** blocksA = ParallelSplit(A, size);
    double*** blocksB = ParallelSplit(B, size1);

    double** D = ParallelMultiply(Sum(blocksA[0], blocksA[3], size / 2), Sum(blocksB[0], blocksB[3], size / 2), size / 2);
    double** H2 = ParallelMultiply(Sum(blocksA[2], blocksA[3], size / 2), blocksB[0], size / 2);
    double** V2 = ParallelMultiply(blocksA[0], Subtract(blocksB[1], blocksB[3], size / 2), size / 2);
    double** V1 = ParallelMultiply(blocksA[3], Subtract(blocksB[2], blocksB[0], size / 2), size / 2);
    double** H1 = ParallelMultiply(Sum(blocksA[0], blocksA[1], size / 2), blocksB[3], size / 2);
    double** D2 = ParallelMultiply(Subtract(blocksA[2], blocksA[0], size / 2), Sum(blocksB[0], blocksB[1], size / 2), size / 2);
    double** D1 = ParallelMultiply(Subtract(blocksA[1], blocksA[3], size / 2), Sum(blocksB[2], blocksB[3], size / 2), size / 2);


    double*** R = new double** [4];
#pragma omp parallel sections
    {
#pragma omp section
        { R[0] = Subtract(Sum(Sum(D, D1, size / 2), V1, size / 2), H1, size / 2); }
#pragma omp section
        { R[1] = Sum(V2, H1, size / 2); }
#pragma omp section
        { R[2] = Sum(V1, H2, size / 2); }
#pragma omp section
        { R[3] = Subtract(Sum(Sum(D, D2, size / 2), V2, size / 2), H2, size / 2); }
    }
    double** result = ParallelCollect(R, size);
    delete[] D, D1, D2, V1, V2, H1, H2, R[0], R[1], R[2], R[3];
    return result;
}


double*** ParallelSplit(double**& matrix, int& size)
{
    if (size % 2 == 1)
    {
        double** temp = new double* [size + 1];
#pragma omp parallel for
        for (int i = 0; i < size + 1; i++)
        {
            temp[i] = new double[size + 1];
            for (int j = 0; j < size + 1; j++)
            {
                if ((i == size) || (j == size))
                {
                    temp[i][j] = 0;
                }
                else
                {
                    temp[i][j] = matrix[i][j];
                }
            }
        }
        matrix = temp;
        size = size + 1;
    }


    double*** blocks = new double** [4];
#pragma omp parallel for
    for (int k = 0; k < 4; k++)
    {
        blocks[k] = new double* [size / 2];
        for (int i = 0; i < size / 2; i++)
        {
            blocks[k][i] = new double[size / 2];
            for (int j = 0; j < size / 2; j++)
            {
                blocks[k][i][j] = matrix[(size / 2) * (k / 2) + i][(size / 2) * (k % 2) + j];
            }
        }
    }

    return blocks;
}


double** ParallelCollect(double*** blocks, int size)
{
    double** matrix = new double* [size];
    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
    }
#pragma omp parallel for
    for (int k = 0; k < 4; k++)
    {
        for (int i = 0; i < size / 2; i++)
        {
            for (int j = 0; j < size / 2; j++)
            {
                matrix[(size / 2) * (k / 2) + i][(size / 2) * (k % 2) + j] = blocks[k][i][j];
            }
        }
    }

    return matrix;
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace DigitsPower
{
    static class AdditionalChainFinder
    {

        static int[] addition_chain_genetic(int[] v, int fit, int popSize, int genNum)
        {
            // v = [7, 31];
            // fit = 8;
            // popSize = 50;
            // genNum = 50000;

            var v_size = size(v);
            var max_number_additional_chain = v[v_size[1] - 1];

            var generation = 0;
            var SelectionProbabilities = rand(popSize, 1);

            var population = InitialPopulation(popSize, max_number_additional_chain);
            var fitness = 1;

            return new[]{1,2,3,4};
        }
        
        static int[,] InitialPopulation(int popSize, int numberGenes)
        {
            return randi(new int[2]{0,1}, popSize, numberGenes);
        }
        
        // TODO evaluate and evaluate_ind functions next 
        
        static int[,] randi(int[] interval, int n, int m)
        {
            Random rand = new Random();
            int [,] array = new int[n,m];
            for (var i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    array[i, j] = rand.Next(interval[0], interval[1]+1);
                }
            }
            return array;
        }
        
        static double[,] rand(int n, int m)
        {
            var matrix = new double[n, m];
            Random random = new Random();
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    matrix[i, j] = random.NextDouble();
                }
            }

            return matrix;
        }

        static int[] size(int[] matrix)
        {
            return new[] {1, matrix.Length};
        }

        static int[] size(int[,] matrix)
        {
            return new[] {matrix.GetLength(0), matrix.GetLength(1)};
        }
    }
}

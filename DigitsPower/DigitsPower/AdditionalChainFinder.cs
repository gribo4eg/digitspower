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
            var fitness = evaluate(population, v);

            var result_check_fitness = check_fitness(fitness, fit);

            while (!result_check_fitness && generation <= genNum)
            {
                int[][][] parents = select(population, SelectionProbabilities);
                population = reproduce(parents, 50, 0.745);
                fitness = evaluate(population, v);
                result_check_fitness = check_fitness(fitness, fit);
                ++generation;
            }
            
            return min(fitness);
        }
        
        static int[][] InitialPopulation (int popSize, int numberGenes)
        {
            return randi(new int[2]{0,1}, popSize, numberGenes);
        }

        static int[][][] select(int[][] pop, double[][] SelectionProbabilities)
        {
            var size_population = size(pop);
            var parents = zeros(2, size_population[0], size_population[1]);

            for (int i = 0; i < size_population[0]; i++)
            {
                var n1 = rand();
                var n2 = rand();
                parents[0][i] = pop[size_population[0]];
                parents[1][i] = pop[size_population[0]];
                for (int j = 0; j < size_population[0]; j++)
                {
                    if (SelectionProbabilities[j][0] >= n1)
                    {
                        parents[0][i] = pop[j];
                    }
                    else if (SelectionProbabilities[j][0] >= n2)
                    {
                        parents[1][i] = pop[j];
                    }
                }
            }

            return parents;
        }

        static int[][] reproduce(int[][][] parents, int mutationDegree, double mutationRate)
        {
            var parents_size = size(parents);
            var result = zeros(parents_size[1], parents_size[2]);
            for (int i = 0; i < parents_size[1]; i++)
            {
                result[i] = crossover_2(parents[0][i], parents[1][i]);
            }

            result = mutate(result, mutationDegree, mutationRate);
            
            return result;
        }

        static int[] crossover_2(int[] parent_1, int[] parent_2)
        {
            var parent_size = size(parent_1);
            var point_1 = randi(new[] {1, parent_size[1] - 1}, 1, 1)[0][0];
            var point_2 = randi(new[] {1, parent_size[1] - 1}, 1, 1)[0][0];

            if (point_1 < point_2)
            {
                return parent_2.Take(point_1)
                    .Concat(parent_1.Skip(point_1).Take(point_2 - point_1 + 1)
                        .Concat(parent_2.Skip(point_2).Take(parent_size[1] - point_2 + 1))).ToArray();
            }
            
            return parent_2.Take(point_2).Concat(parent_1.Skip(point_2).Take(point_1 - point_2 + 1)
                    .Concat(parent_2.Skip(point_1).Take(parent_size[1] - point_1 + 1))).ToArray();
            
        }

        static int[][] mutate(int[][] pop, int mutationDegree, double mutationRate)
        {
            var size_population = size(pop);
            if (Convert.ToBoolean(mutationRate) && Convert.ToBoolean(mutationDegree))
            {
                for (int i = 0; i < size_population[0]; i++)
                {
                    double n = rand();
                    if (n <= mutationRate)
                    {
                        for (int j = 0; j < mutationDegree; j++)
                        {
                            var gene = randi(new[] {1, size_population[1]}, 1, 1)[0][0];
                            pop[i][gene] = (pop[i][gene] + 1) % 2;
                        }
                    }
                }
            }

            return pop;
        }

        static int[] evaluate (int[][] population, int[] v)
        {
            var size_population = size(population);
            //var fitness = zeros(1, size_population[0]);
            var fitness = new int[size_population[0]];

            for (int i = 0; i < size_population[0]; i++)
            {
                fitness[i] = evaluate_ind(population[i], v, size_population[1]);
            }
            
            return fitness;
        }

        static int evaluate_ind(int[] s, int[] v, int n)
        {
            int largePenalty = int.MaxValue;
            int fitness = 0;

            for (int i = 0; i < n; i++)
            {
                if (Convert.ToBoolean(s[i]))
                {
                    ++fitness;
                }

                if (!isInclude(s, v) || (Convert.ToBoolean(s[i]) && !isSumOfPrevious(s, i)))
                {
                    fitness += largePenalty;
                    break;
                }
            }
            return fitness;
        }

        static bool check_fitness(int[] fitness, int value)
        {
            bool out_v = false;
            var fitness_size = size(fitness);
            for (int i = 0; i < fitness_size[1]; i++)
            {
                if (fitness[i] <= value)
                {
                    out_v = true;
                    break;
                }
            }

            return out_v;
        }

        static bool isInclude(int[] s, int[] v)
        {
            bool out_v = true;
            var size_v = size(v);
            for (int i = 0; i < size_v[1]; i++)
            {
                if (Convert.ToBoolean(s[v[i]])) continue;
                out_v = false;
                break;
            }

            return out_v;
        }

        static bool isSumOfPrevious(int[] s, int i)
        {
            bool out_v = false;
            for (int j = 0; j < i-1; j++)
            {
                for (int k = 0; k < i-1; k++)
                {
                    if ((j + k == 1) && (s[j] == 1) && (s[k] == 1))
                    {
                        return true;
                    }
                }
            }

            return out_v;
        }

        static int[] min(int[] array)
        {
            int min = array.Min();
            int index = Array.IndexOf(array, min);
            return new[] {min, index};
        }
        
        static int[][] randi (int[] interval, int n, int m)
        {
            Random rand = new Random();
            int [][] array = new int[n][];
            for (var i = 0; i < n; i++)
            {
                array[i] = new int[m];
                for (int j = 0; j < m; j++)
                {
                    array[i][j] = rand.Next(interval[0], interval[1]+1);
                }
            }
            return array;
        }
        
        static double[][] rand (int n, int m)
        {
            double[][] matrix = new double[n][];
            Random random = new Random();
            for (int i = 0; i < n; i++)
            {
                matrix[i] = new double[m];
                for (int j = 0; j < m; j++)
                {
                    matrix[i][j] = random.NextDouble();
                }
            }

            return matrix;
        }

        static double rand()
        {
            Random rand = new Random();

            return rand.NextDouble();
        }

        static int[] size (int[] matrix)
        {
            return new[] {1, matrix.Length};
        }

        static int[] size (int[][] matrix)
        {
            return new[] {matrix.Length, matrix[0].Length};
        }

        static int[] size(int[][][] matrix)
        {
            return new[] {matrix.Length, matrix[0].Length, matrix[0][0].Length};
        }

        static int[][] zeros(int n, int m)
        {
            int[][] matrix = new int[n][];
            for (int i = 0; i < n; i++)
            {
                matrix[i] = new int[m];
                for (int j = 0; j < m; j++)
                {
                    matrix[i][j] = 0;
                }
            }

            return matrix;
        }

        static int[][][] zeros(int n, int m, int l)
        {
            int[][][] matrix = new int[n][][];
            for (int i = 0; i < n; i++)
            {
                matrix[i] = new int[m][];
                for (int j = 0; j < m; j++)
                {
                    matrix[i][j] = new int[l];
                    for (int k = 0; k < l; k++)
                    {
                        matrix[i][j][k] = 0;
                    }
                }
            }

            return matrix;
        }
    }
}

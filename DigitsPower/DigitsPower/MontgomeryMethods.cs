using System;
using System.Numerics;
using static System.Math;
using static DigitsPower.HelpMethods;



namespace DigitsPower
{
    public static class MontgomeryMethods
    {

        public static MyList<BigInteger> toMontgomeryDomain(ref BigInteger a, ref BigInteger b, BigInteger N)
        {
            String n_str = ConvToBinary(N);
            int m = n_str.Length;
            BigInteger r = BigInteger.One << m;
            BigInteger r_inv = Euclid_2_1(N, r);
            BigInteger n_shtrih = (r * r_inv - 1) / N;
            BigInteger b_module = r - 1; // маска для виконання операцій за модулем b

            MyList<BigInteger> result = new MyList<BigInteger>();
            result.Add(n_shtrih);
            result.Add(b_module);
            result.Add(m);

            /*
            // обчислення a_Montgomery (a = a * r % N)
            a = a * r % N;
            // обчислення b_Montgomery (b = b * r % N)
            b = b * r % N;
            */

            BigInteger r_sqr = r * r % N;
            result.Add(r_sqr * r % N);

            // обчислення a_Montgomery (a = a * r % N)
            a = MontgomeryMultDomain(a, r_sqr, N, result); // це буде a*(r^2), потім отримуємо a*(r^2)*(r^-1) = a * r mod N

            // обчислення b_Montgomery (b = b * r % N)
            b = MontgomeryMultDomain(b, r_sqr, N, result);

            return result;
        }

        public static BigInteger outMontgomeryDomain(BigInteger result, BigInteger N, MyList<BigInteger> parameters)
        {
            return MontgomeryReduction(result, N, parameters);
        }

        public static BigInteger MontgomeryReduction(BigInteger result, BigInteger N, MyList<BigInteger> parameters)
        {
            result = (result + ((result * parameters[0]) & parameters[1]) * N) >> (int)parameters[2];
            if (result >= N)
                result -= N;

            return result;
        }
        
        public static BigInteger MontgomeryMultDomain(BigInteger a, BigInteger b, BigInteger N, MyList<BigInteger> parameters)
        {
            return MontgomeryReduction(a * b, N, parameters);
        }

        public static BigInteger MontgomeryInverse(BigInteger mod, BigInteger found, MyList<BigInteger> parameters)
        {
            BigInteger inv = Euclid_2_1(mod, found);
            return MontgomeryMultDomain(inv, parameters[3], mod, parameters); // inv домножити на r^3 використовуючи множення Монтгомері
        }

        /*
        public static BigInteger MontgomeryInverse(BigInteger mod, BigInteger found)
        {
            // насправді це метод 2.4

            BigInteger u, v, B, D, y, t1, t2, q, d, inv;
            u = mod;
            v = found;

            B = 0;
            D = 1;
            int k = 0;
            while (v != 0)
            {
                if (u % 2 == 0)
                {
                    u = u >> 1;
                    D = D << 1;
                }
                else if (v % 2 == 0)
                {
                    v = v >> 1;
                    B = B << 1;
                }
                else if (u > v)
                {
                    u = (u - v) >> 1;
                    B = B + D;
                    D = D << 1;
                }
                else
                {
                    v = (v - u) >> 1;
                    D = D + B;
                    B = B << 1;
                }

                k++;
               
            }

            //d = u; y = B;

            if (B >= mod) B = B - mod;

            return mod - B;
        }
*/
 
    }
}
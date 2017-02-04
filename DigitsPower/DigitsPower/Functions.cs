using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;
using System.Diagnostics;
using System.Security.Cryptography;

namespace MultInverse
{
    public class Functions
    {
        public static int rand(int l, int h)
        {
            Random r = new Random();
            int x = r.Next(l, h);
            return x;
        }

        public static BigInteger random_max(int n)
        {
            var rng = new RNGCryptoServiceProvider();
            byte[] bytes = new byte[n + 1];
            rng.GetBytes(bytes);
            bytes[n - 1] |= 0x80;
            bytes[n] = 0;
            bytes[0] |= 1;
            BigInteger b = new BigInteger(bytes);

            return b;
        }

        public static string ToBin(BigInteger dec)
        {
            string BinResult = "";
            string a;

            while (dec > 0)
            {
                BinResult += dec % 2;
                dec >>= 1;
            }
            return a = new string(BinResult.Reverse().ToArray());
        }

        public static List<BigInteger> NAFw(BigInteger k, int w)
        {
            List<BigInteger> mas_k = new List<BigInteger>();
            int i = 0;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    BigInteger temp = (k % (int)Math.Pow(2, w));
                    if (temp >= BigInteger.Pow(2, w - 1))
                        mas_k.Add(temp - (int)Math.Pow(2, w));
                    else
                        mas_k.Add(temp);

                    k = k - mas_k[i];
                }
                else
                    mas_k.Add(0);

                k = k / 2;
                i++;
            }
            return mas_k;

        }

        public static List<BigInteger> NAF(BigInteger k)
        {
            List<BigInteger> mas_k = new List<BigInteger>();
            int i = 0;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    mas_k.Add(2 - (k % 4));
                    k = k - mas_k[i];
                }
                else
                    mas_k.Add(0);


                k = k / 2;
                i++;
            }
            return mas_k;
        }

        public static void Find_the_largest_t(List<BigInteger> mas_k, int i, int w, out int max, out int max_j)
        {

            int j = i;
            BigInteger pow = 1;
            BigInteger temp = pow * mas_k[i];
            max = (int)temp;
            max_j = i;
            while (j - i + 1 < w && j < mas_k.Count - 1)
            {
                j++;
                pow = pow * 2;
                temp = temp + pow * mas_k[j];
                if (temp % 2 != 0)
                {
                    max_j = j;
                    max = (int)temp;
                }
            }

            max_j = max_j - i + 1;

        }

        public static void Find_the_largest_t_10(List<BigInteger> mas_k, int i, int w, out int max, out int max_j)
        {

            int j = i;
            int pow;
            BigInteger temp = mas_k[i];
            max = (int)temp;
            max_j = i;
            while (i - j < w - 1 && j > 0)
            {
                pow = 1;
                temp = 0;
                j--;
                for (int t = j; t <= i; t++)
                {
                    temp = temp + pow * mas_k[t];
                    pow = pow * 2;
                }

                if (temp % 2 != 0)
                {
                    max_j = j;
                    max = (int)temp;
                }
            }

            max_j = i - max_j + 1;

        }

        public static void Extended_Euclid(BigInteger a, BigInteger b, out BigInteger d, out BigInteger inv)
        {
            BigInteger u, v, A, C, t1, t2, y;

            if (a < b)
                b = b % a;

            u = a;
            v = b;
            A = 1;
            C = 0;

            while (v != 0)
            {

                BigInteger q = u / v;
                if (q * v > u) q--;
                t1 = u - q * v;
                t2 = A - q * C;

                u = v;
                A = C;

                v = t1;
                C = t2;

            }

            //x = A;
            if (b == 0) { d = u; inv = 0; }
            else
            {
                y = (u - A * a) / b;
                d = u;


                if (y >= 0)
                    inv = y;
                else
                    inv = y + a;
            }

        }

        public static BigInteger JACOBI(BigInteger a, BigInteger n)
        {
            int s = 0;

            BigInteger a1, b = a, e = 0, m, n1;


            if (a == 0) return 0;

            if (a == 1) return 1;

            while ((b % 2) == 0)
            {

                b >>= 1;
                e++;
            }

            a1 = b;

            m = n % 8;

            if ((e % 2) == 0) s = 1;

            else if (m == 1 || m == 7) s = +1;

            else if (m == 3 || m == 5) s = -1;

            if (n % 4 == 3 && a1 % 4 == 3) s = (-1) * s;

            if (a1 != 1) n1 = n % a1; else n1 = 1;

            return s * JACOBI(n1, a1);

        }

        public static BigInteger square_root_mod(BigInteger a, BigInteger p)

/* returns the square root of a modulo an odd prime p

   if it exists 0 otherwise */
        {

            BigInteger ai, b = 0, c, d, e, i, r, s = 0, t = p - 1;


            /* is a quadratic nonresidue */

            if (JACOBI(a, p) == -1) return 0;

            /* find quadratic nonresidue */

            do
            {

                do
                {
                    b = rand(1, 32000) % p;
                }
                while (b == 0);
            }
            while (JACOBI(b, p) != -1);

            /* write p - 1 = 2 ^ s * t for odd t */

            while ((t % 2) == 0)
            {
                s++;
                t >>= 1;
            }

            BigInteger d_, inv;

            double time;
            Lehmer_.Lehmer_long(p, a, out d_, out inv, out time);
            ai = inv;

            c = BigInteger.ModPow(b, t, p);

            r = BigInteger.ModPow(a, (t + 1) / 2, p);

            for (i = 1; i < s; i++)
            {

                e = BigInteger.ModPow(2, s - i - 1, p);

                d = BigInteger.ModPow((((r * r) % p) * ai) % p, e, p);

                if (d == p - 1) r = r * c % p;

                c = c * c % p;

            }

            return r;

        }

        public static BigInteger GreatestPrimeDivisor(BigInteger x)
        {
            BigInteger res = 0, temp = 0;
            bool flag = true;

            for (BigInteger i = x << 1; i > 2; i--)
            {
                if (x % i == 0)
                {
                    for (BigInteger j = 2; j < i; j++)
                    {
                        if (i % j == 0)
                        {
                            flag = false;
                            break;
                        }
                        else flag = true;
                    }
                    if (flag == true) res = i;
                }
                if (res != 0) break;
            }
            if (res == 0) res = 2; //0 polu4aetsya, esli x = 2,4,8...
            return res;
        }

        public static BigInteger sqrt_int(BigInteger b)
        {
            BigInteger x = 1;
            bool decreased = false;
            for (; ; )
            {
                BigInteger nx = (x + b / x) >> 1;
                if (x == nx || nx > x && decreased) break;
                decreased = nx < x;
                x = nx;
            }
            return x;
        }

        public static BigInteger Pow(BigInteger value, BigInteger exponent)
        {
            BigInteger _base = 1;
            for (BigInteger i = 0; i < exponent; i++)
            {
                _base *= value;
            }
            return _base;
        }

        public static double log_dif_base(BigInteger _base, double argument)
        {
            return (Math.Log(argument) / BigInteger.Log(_base));
        }

        public static BigInteger get_number_bit(double value, BigInteger _base)
        {
            BigInteger number_bit = (BigInteger)Math.Floor(log_dif_base(_base, (double)value));
            if ((log_dif_base(_base, (double)value) > (double)number_bit) || number_bit == 0)
                number_bit += 1;
            return number_bit;
        }

        public static void Re_Compute_a_b(ref BigInteger a, ref BigInteger b, BigInteger a_max, BigInteger b_max, double k)
        {
            if (a > a_max)
            {
                BigInteger temp = get_number_bit((double)Pow(2, a - a_max), 3) - 1;
                b = b + temp;

                if (b > b_max)
                {
                    b = b_max;
                }
                if (a_max > 0)
                {
                    temp = a_max - 1;
                }
                else temp = 0;
                if ((Math.Abs(k - (double)Pow(2, a_max) * (double)Pow(3, b)) < Math.Abs(k - (double)Pow(2, temp) * (double)Pow(3, b + 1))) || b == b_max)
                {
                    a = a_max;
                }
                else
                {
                    a = temp;
                    b = b + 1;
                }
            }
        }

        public static void Best_Approximation_1(double k, BigInteger a_max, BigInteger b_max, out BigInteger a, out BigInteger b)
        {
            BigInteger min_x = a_max;
            double tmp = (double)log_dif_base(3, 2);
            BigInteger y = (BigInteger)Math.Round(((double)(-min_x)) * (double)log_dif_base(3, 2) + (double)log_dif_base(3, k));
            if (y > b_max)
            {
                y = b_max;
            }
            double min_delta = Math.Abs(((double)y + ((double)min_x * log_dif_base(3, 2)) - log_dif_base(3, k)));
            for (BigInteger x = 1; x < a_max; x++)
            {
                y = (BigInteger)Math.Round(((double)(-x)) * (double)log_dif_base(3, 2) + (double)log_dif_base(3, k));
                if (y > b_max)
                {
                    y = b_max;
                }
                double delta = Math.Abs(((double)y + (double)x * log_dif_base(3, 2) - log_dif_base(3, k)));
                if (min_delta > delta)
                {
                    min_x = x;
                    min_delta = delta;
                }
            }

            a = min_x;
            b = (BigInteger)Math.Round(((double)(-min_x)) * (double)log_dif_base(3, 2) + (double)log_dif_base(3, k));
            if (b > b_max)
            {
                b = b_max;
            }
        }

        public static void Best_Approximation_2(double k, BigInteger a_max, BigInteger b_max, out BigInteger a, out BigInteger b)
        {
            a = 0;
            b = 0;
            BigInteger legth_k = get_number_bit(k, 2);
            BigInteger[,] PreComputation = new BigInteger[(int)b_max + 1, 3];
            int i;
            for (i = 0; i <= b_max; i++)
            {
                PreComputation[i, 0] = i;
                PreComputation[i, 1] = BigInteger.Pow(3, i);
                BigInteger temp = get_number_bit((double)PreComputation[i, 1], 2);
                PreComputation[i, 2] = (PreComputation[i, 1] << (int)(legth_k - temp));
            }

            for (i = 0; i <= b_max; i++)
            {
                int i_min = i;
                BigInteger min = PreComputation[i, 2];
                for (int u = i + 1; u <= b_max; u++)
                {
                    if (min > PreComputation[u, 2])
                    {
                        i_min = u;
                        min = PreComputation[u, 2];
                    }
                }
                for (int j = 0; j < 3; j++)
                {
                    BigInteger temp = PreComputation[i_min, j];
                    PreComputation[i_min, j] = PreComputation[i, j];
                    PreComputation[i, j] = temp;
                }
            }
            i = (int)b_max + 1;
            int length_1 = 0;
            int length_max = 0;

            while (i > 0 && length_1 >= length_max)
            {
                int j = 0;
                length_max = length_1;
                string str1 = ToBin((BigInteger)k);
                string str2 = ToBin(PreComputation[i - 1, 2]);
                while (j < (int)legth_k - 1 && (str2[j] == str1[j]))
                {
                    j = j + 1;
                }
                if (j != 0)
                {
                    length_1 = (int)legth_k - ((int)legth_k - j);
                }
                else length_1 = (int)legth_k;
                i = i - 1;
            }
            if (length_1 < length_max)
            {
                i = i + 2;
            }
            else i = 1;

            BigInteger b1 = PreComputation[i - 1, 0];
            BigInteger a1 = legth_k - get_number_bit((double)PreComputation[i - 1, 1], 2);
            if (a1 < 0)
            {
                a1 = 0;
            }
            BigInteger a2, b2;
            if (i < b_max + 1)
            {
                b2 = PreComputation[i, 0];
                a2 = legth_k - get_number_bit((double)PreComputation[i, 1], 2);
                if (a2 < 0)
                {
                    a2 = 0;
                }
            }
            else
            {
                b2 = 0;
                a2 = 0;
            }

            Re_Compute_a_b(ref a1, ref b1, a_max, b_max, k);
            Re_Compute_a_b(ref a2, ref b2, a_max, b_max, k);
            if (Math.Abs(k - (double)Pow(2, a1) * (double)Pow(3, b1)) < Math.Abs(k - (double)Pow(2, a2) * (double)Pow(3, b2)))
            {
                a = a1;
                b = b1;
            }
            else
            {
                a = a2;
                b = b2;
            }
            if ((a != a_max) && (Math.Abs(k - (double)Pow(2, a + 1) * (double)Pow(3, b)) < Math.Abs(k - (double)Pow(2, a) * (double)Pow(3, b))))
            {
                a = a + 1;
            }
        }
        
        public static BigInteger[,] Convert_to_DBNS_1(double k, BigInteger a_max, BigInteger b_max)
        {
            List<BigInteger[]> mass_k_l = new List<BigInteger[]>();
            int i = 0;
            int s = 1;
            BigInteger a, b;
            while (k > 0)
            {
                i++;
                Best_Approximation_1(k, a_max, b_max, out a, out b);
                a_max = a;
                b_max = b;
                double z = Math.Pow(2, (double)a) * Math.Pow(3, (double)b);
                mass_k_l.Add(new BigInteger[3]{s, a, b});
                if (k < z)
                {
                    s = -s;
                }
                k = Math.Abs(k - z);   
            }
            BigInteger[,] mass_k = new BigInteger[mass_k_l.Count, 3];
            for (int j = 0; j < mass_k_l.Count; j++)
            {
                mass_k[j,0] = mass_k_l[j][0];
                mass_k[j, 1] = mass_k_l[j][1];
                mass_k[j, 2] = mass_k_l[j][2];
            }
            return mass_k;
        }
        
        public static BigInteger[,] Convert_to_DBNS_2(double k, BigInteger a_max, BigInteger b_max)
        {
            List<BigInteger[]> mass_k_l = new List<BigInteger[]>();
            int i = 0;
            int s = 1;
            BigInteger a, b;
            while (k > 0)
            {
                i++;
                Best_Approximation_2(k, a_max, b_max, out a, out b);
                a_max = a;
                b_max = b;
                double z = Math.Pow(2, (double)a) * Math.Pow(3, (double)b);
                mass_k_l.Add(new BigInteger[3] { s, a, b });
                if (k < z)
                {
                    s = -s;
                }
                k = Math.Abs(k - z);
            }
            BigInteger[,] mass_k = new BigInteger[mass_k_l.Count, 3];
            for (int j = 0; j < mass_k_l.Count; j++)
            {
                mass_k[j, 0] = mass_k_l[j][0];
                mass_k[j, 1] = mass_k_l[j][1];
                mass_k[j, 2] = mass_k_l[j][2];
            }
            return mass_k;
        }
    }
}

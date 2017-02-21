using System;
using System.Collections.Generic;
using System.Numerics;
using static DigitsPower.PowFunctions;
using System.Linq;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using static DigitsPower.HelpMethods;
namespace DigitsPower
{
    static class HelpMethods
    {
        #region MultiplyMethods
        public static BigInteger SimpleMultiply(BigInteger a, BigInteger b, BigInteger m)
        {
            return (a * b) % m;
        }

        #endregion
        #region BinaryMethods
        //public static BigInteger Point_Multiplication_Affine_Coord_21(BigInteger found, BigInteger pow, BigInteger mod)
        //{
        //    List<int[]> mas_k = PowFunctions.ToDBNS2RL(pow);

        //    BigInteger res = 1;
        //    BigInteger t = found;



        //    for (int i = 0; i < mas_k.Count; i++)
        //    {
        //        for (int j = 0; j < mas_k[i][1]; j++)
        //            t = (t * t) % mod;
        //        for (int j = 0; j < mas_k[i][2]; j++)
        //            t = ((t * t) % mod * t) % mod;

        //        if (mas_k[i][0] == 1)
        //            res = res * t % mod;
        //        else if (mas_k[i][0] == -1)
        //            res = res * inv(t, mod) % mod;
        //    }
        //    return res;
        //}

        public static List<BigInteger[]> Convert_to_DBNS_1_1(BigInteger k, int a_max, int b_max)
        {
            List<BigInteger[]> mass_k_l = new List<BigInteger[]>();
            int s = 1;
            List<int> app;
            int a = a_max;
            int b = b_max;
            while (k > 0)
            {
                app = Best_Approximation_1_1(k, a, b);
                a = app[0];
                b = app[1];
                BigInteger z = (BigInteger.One << a) * BigInteger.Pow(3, b);
                mass_k_l.Add(new BigInteger[3] { s, a, b });
                if (k < z)
                    s = -s;

                k = BigInteger.Abs(k - z);
            }

            return mass_k_l;
        }

        public static List<BigInteger[]> Convert_to_DBNS_1_2(BigInteger k, int a_max, ref int b_max)
        {
            List<BigInteger[]> mass_k_l = new List<BigInteger[]>();
            BigInteger i = 0;
            BigInteger s = 1;
            List<int> app;

            BigInteger temp = 1;
            BigInteger three = 3;
            SortedList<string, BigInteger> PreComputation = new SortedList<string, BigInteger>();

            double log_3 = BigInteger.Log(k, 3);

            if (b_max > log_3)
                b_max = (int)Math.Floor(log_3);

            for (i = 0; i <= b_max; i++)
            {
                PreComputation.Add(ConvToBinary(temp), i);
                temp *= three;
            }

            int a = a_max;
            int b = b_max, b_temp = b, index;
            while (k > 0)
            {
                log_3 = BigInteger.Log(k, 3);
                if (b > log_3)
                    b = (int)Math.Round(log_3);

                for (i = b_temp; i > b; i--)
                {
                    index = PreComputation.IndexOfValue(i);
                    PreComputation.RemoveAt(index);
                }

                b_temp = b;
                app = Best_Approximation_1_2(k, PreComputation, a, b);

                a = app[0];
                b = app[1];
                BigInteger z = (BigInteger.One << a) * BigInteger.Pow(3, b);
                mass_k_l.Add(new BigInteger[3] { s, a, b });
                if (k < z)
                    s = -s;

                k = BigInteger.Abs(k - z);
            }

            return mass_k_l;
        }

        public static List<int> Best_Approximation_1_1(BigInteger k, int a_max, int b_max)
        {
            int a;
            int b;
            int min_x = a_max;
            double log_3_k = log_dif_base(3, k);
            double log_3_2 = log_dif_base(3, 2);
            double temp = -min_x * log_3_2 + log_3_k;
            double min_temp = temp;
            int y = (int)Math.Round(temp);
            if (y > b_max)
                y = b_max;
            else if (y < 0)
                y = 0;

            double min_delta = Math.Abs(y - temp);
            for (int x = 0; x < a_max; x++)
            {
                temp = -x * log_3_2 + log_3_k;
                y = (int)Math.Round(temp);
                if (y > b_max)
                    y = b_max;
                else if (y < 0)
                    y = 0;

                double delta = Math.Abs(y - temp);
                if (min_delta > delta)
                {
                    min_x = x;
                    min_temp = temp;
                    min_delta = delta;
                }
            }

            a = min_x;
            b = (int)Math.Round(min_temp);
            if (b > b_max)
                b = b_max;

            List<int> r = new List<int>();
            r.Add(a);
            r.Add(b);
            return r;
        }

        public static List<int> Best_Approximation_1_2(BigInteger k, SortedList<string, BigInteger> PreComputation, int a_max, int b_max)
        {
            int a = 0;
            int b = 0;

            if (k > 1)
            {
                string temp_str = ConvToBinary(k);

                bool fl = false;
                int index_k;
                if (!PreComputation.ContainsKey(temp_str))
                {
                    PreComputation.Add(temp_str, k);
                    fl = true;
                }
                index_k = PreComputation.IndexOfKey(temp_str);

                int a1, a2;
                int b1, b2;

                if (index_k > 0)
                {
                    b1 = (int)PreComputation.ElementAt(index_k - 1).Value;
                    a1 = temp_str.Length - PreComputation.ElementAt(index_k - 1).Key.Length;
                    a1 = change_value(a_max, a1);
                }
                else
                {
                    b1 = 0;
                    a1 = change_value(a_max, temp_str.Length);
                }

                if (index_k < PreComputation.Count - 1)
                {
                    b2 = (int)PreComputation.ElementAt(index_k + 1).Value;
                    a2 = temp_str.Length - PreComputation.ElementAt(index_k + 1).Key.Length;
                    a2 = change_value(a_max, a2);
                }
                else
                {
                    b2 = 0;
                    a2 = change_value(a_max, temp_str.Length);
                }


                int a_dif, b_dif;
                BigInteger power_a1, power_a2, power_b1, power_b2;
                if (a1 > a2)
                {
                    a_dif = a1 - a2;
                    power_a2 = BigInteger.One << a2;
                    power_a1 = power_a2 << a_dif;
                }
                else
                {
                    a_dif = a2 - a1;
                    power_a1 = BigInteger.One << a1;
                    power_a2 = power_a1 << a_dif;
                }

                if (b1 > b2)
                {
                    b_dif = b1 - b2;
                    power_b2 = BigInteger.Pow(3, b2);
                    power_b1 = power_b2 << b_dif;
                }
                else
                {
                    b_dif = b2 - b1;
                    power_b1 = BigInteger.Pow(3, b1);
                    power_b2 = power_b1 << b_dif;
                }

                BigInteger delta_1 = BigInteger.Abs(k - power_a1 * power_b1);
                BigInteger delta_2 = BigInteger.Abs(k - power_a2 * power_b2);

                if (delta_1 < delta_2)
                {
                    a = a1;
                    b = b1;
                }
                else
                {
                    a = a2;
                    b = b2;
                }

                if (fl)
                    PreComputation.RemoveAt(index_k);
            }

            List<int> r = new List<int>();
            r.Add(a);
            r.Add(b);

            return r;


            //int next = 0;

            /*
                long legth_k = get_number_bit(k, 2);
            long[,] PreComputation = new long[b_max + 1, 3];
            //long i;
            for (i = 0; i <= b_max; i++)
            {
                PreComputation[i, 0] = i;
                PreComputation[i, 1] = (long)Math.Pow(3, i);
                long temp = get_number_bit(PreComputation[i, 1], 2);
                if (temp > legth_k)
                {
                    b_max = i - 1;
                    break;
                }
                PreComputation[i, 2] = (PreComputation[i, 1] << (int)(legth_k - temp));
            }

            for (i = 0; i <= b_max; i++)
            {
                long i_min = i;
                long min = PreComputation[i, 2];
                for (long u = i + 1; u <= b_max; u++)
                {
                    if (min > PreComputation[u, 2])
                    {
                        i_min = u;
                        min = PreComputation[u, 2];
                    }
                }
                for (long j = 0; j < 3; j++)
                {
                    long temp = PreComputation[i_min, j];
                    PreComputation[i_min, j] = PreComputation[i, j];
                    PreComputation[i, j] = temp;
                }
            }
            i = b_max + 1;
            long length_1 = 0;
            long length_max = 0;

            while (i > 0 && length_1 >= length_max)
            {
                int j = 0;
                length_max = length_1;
                string str1 = bif.ToBin(k);
                string str2 = Convert.ToString((PreComputation[i - 1, 2]), 2);
                while (j < legth_k - 1 && j < Math.Min(str2.Length, str1.Length) && (str2[j] == str1[j]))
                {
                    j = j + 1;
                }
                if (j != 0)
                {
                    length_1 = legth_k - (legth_k - j);
                }
                else length_1 = legth_k;
                i = i - 1;
            }
            if (length_1 < length_max)
            {
                i = i + 2;
            }
            else i = 1;

            long b1 = PreComputation[i - 1, 0];
            long a1 = legth_k - get_number_bit(PreComputation[i - 1, 1], 2);
            if (a1 < 0)
            {
                a1 = 0;
            }
            long a2, b2;
            if (i < b_max + 1)
            {
                b2 = PreComputation[i, 0];
                a2 = legth_k - get_number_bit(PreComputation[i, 1], 2);
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
            if (Math.Abs((double)k - Math.Pow(2, a1) * Math.Pow(3, b1)) < Math.Abs((double)k - Math.Pow(2, a2) * Math.Pow(3, b2)))
            {
                a = a1;
                b = b1;
            }
            else
            {
                a = a2;
                b = b2;
            }
            if ((a != a_max) && (Math.Abs((double)k - Math.Pow(2, a + 1) * Math.Pow(3, b)) < Math.Abs((double)k - Math.Pow(2, a) * Math.Pow(3, b))))
            {
                a = a + 1;
            }


            List<int> r = new List<int>();
            r.Add(a);
            r.Add(b);

            return r;*/
        }

        public static int change_value(int max, int value)
        {
            if (value <= max)
                return value;
            else
                return max;
        }

        public static double log_dif_base(BigInteger _base, BigInteger argument)
        {
            return (BigInteger.Log(argument) / BigInteger.Log(_base));
        }

        public static long get_number_bit(BigInteger value, long _base)
        {
            long number_bit = (long)Math.Floor(bif.log_dif_base(value, _base));
            if (((long)bif.log_dif_base(value, _base) > number_bit) || number_bit == 0)
                number_bit++;
            return number_bit;
        }

        public static void Re_Compute_a_b(ref long a, ref long b, long a_max, long b_max, BigInteger k)
        {
            if (a > a_max)
            {
                long temp = get_number_bit((long)Math.Pow(2, a - a_max), 3) - 1;
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
                if ((Math.Abs((double)k - Math.Pow(2, a_max) * Math.Pow(3, b)) < Math.Abs((double)k - Math.Pow(2, temp) * Math.Pow(3, (b + 1)))) || b == b_max)
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

        public static BigInteger Euclid_2_1(BigInteger mod, BigInteger found)
        {
            // насправді це метод 2.4

            BigInteger u, v, B, D, y, t1, t2, q, d, inv;
            u = mod;
            v = found;

            B = 0;
            D = 1;
            while (v != 0)
            {
                q = u / v;
                t1 = u % v;
                t2 = B - q * D;
                u = v;
                B = D;

                v = t1;
                D = t2;
            }
            d = u; y = B;

            if (y >= 0) inv = y;
            else inv = y + mod;

            /*
            BigInteger u, v, A, B, C, D, y, t1, t2, t3, q, d, inv;

            u = mod;
            v = found;

            A = 1;
            B = 0;
            C = 0;
            D = 1;

            while (v != 0)
            {
                q = u / v;
                t1 = u - q * v;
                t2 = A - q * C;
                t3 = B - q * D;

                u = v;
                A = C;
                B = D;

                v = t1;
                C = t2;
                D = t3;
            }
            d = u; y = B;

            if (y >= 0) inv = y;
            else inv = y + mod;
            */

            return inv;

        }

        public static BigInteger NSD(BigInteger mod, BigInteger found)
        {
            BigInteger u, v, A, B, C, D, t1, t2, t3, q, d;

            u = mod;
            v = found;

            A = 1;
            B = 0;
            C = 0;
            D = 1;

            while (v != 0)
            {
                q = u / v;
                t1 = u - q * v;
                t2 = A - q * C;
                t3 = B - q * D;

                u = v;
                A = C;
                B = D;

                v = t1;
                C = t2;
                D = t3;
            }
            d = u;

            return d;
        }

        public static BigInteger TwoPow(BigInteger pow)
        {
            BigInteger result = 1;
            for (BigInteger i = 0; i < pow; i++)
            {
                result = result << 1;
            }
            return result;
        }
        public static string ConvToBinary(BigInteger dec)
        {
            string BinResult = "";

            while (dec > 0)
            {
                BinResult = dec % 2 + BinResult;
                dec >>= 1;
            }

            return BinResult;
        }
        public static MyList<int> ToNAF(BigInteger pow)
        {
            MyList<int> mas_k = new MyList<int>();
            BigInteger k = pow;
            Int32 temp_1, temp_2;
            int i = 0;
            while (k >= 1)
            {
                temp_1 = (Int32)(k & 3); // k & 3 це те саме, що k % 4
                temp_2 = temp_1 & 1;     // k & 1 це те саме, що k % 2
                if (temp_2 != 0)
                {
                    mas_k.Add((int)(2 - temp_1));
                    k = k - mas_k[i];
                }
                else
                    mas_k.Add(0);

                k = k >> 1;
                i++;
            }
            return mas_k;

            /*
            MyList<int> res = new MyList<int>();
            BigInteger k = pow;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    res.Insert(0, 2 - (int)(k % 4));
                    k = k - res[0];
                }
                else
                    res.Insert(0, 0);
                k = k >> 1;
            }
            return res;
            */
        }
        public static MyList<int[]> ToDBNS2RL(BigInteger pow)
        {
            var x = new MyList<int[]>();
            BigInteger v = pow;
            int i = 0;
            while (v != 0)
            {
                x.Insert(i, new int[3]);
                int a = 0;
                while (v % 2 == 0 && v != 0)
                {
                    a++;
                    v >>= 1;
                }
                int b = 0;
                while (v % 3 == 0 && v != 0)
                {
                    b++;
                    v /= 3;
                }
                x[i][1] = a;
                x[i][2] = b;

                v--;
                if (v != 0)
                {
                    if (v % 6 != 0)
                    {
                        v += 2;
                        x[i][0] = -1;
                    }
                    else
                        x[i][0] = 1;
                }
                else
                    x[i][0] = 1;
                i++;
            }
            return x;
        }

        public static MyList<int[]> ToDBNS2LR(BigInteger pow)
        {
            MyList<int[]> x = new MyList<int[]>();
            x.Add(new int[3]);
            int i = 0;
            BigInteger v = pow;
            x[i][0] = 0;
            while (v != 0)
            {
                int a = 0;
                while (v % 2 == 0 && v != 0)
                {
                    a++;
                    v >>= 1;
                }

                int b = 0;
                while (v % 3 == 0 && v != 0)
                {
                    b++;
                    v = v / 3;
                }
                x[i][1] = a;
                x[i][2] = b;
                i++;
                x.Add(new int[3]);
                v--;
                if (v != 0)
                {
                    if (v % 6 != 0)
                    {
                        v = v + 2;
                        x[i][0] = -1;
                    }
                    else
                        x[i][0] = 1;
                }
            }
            return x;
        }

        #endregion


        #region WindowsMethods

        public static MyList<int> FindLargest1(MyList<int> x, int i, int w)
        {
            // виділяємо вікно в масиві х від молодшого індексу і до старших індексів
            // використовуємо коли рухаємось від 0 до x.Count - 1

            int j = i;
            int pow = 1;
            int temp = x[i];

            int max_j = i;
            int max = x[i];
            while (j - i + 1 < w && j < x.Count - 1)
            {
                j++;

                pow = pow << 1;
                temp = temp + pow * x[j];

                if (x[j] != 0)
                {
                    max_j = j;
                    max = temp;
                }
            }

            max_j = max_j - i + 1;

            MyList<int> r = new MyList<int>();
            r.Add(max);
            r.Add(max_j);
            return r;

        }
        public static MyList<int> FindLargest2(MyList<int> x, int i, int w)
        {
            // виділяємо вікно в масиві х від старшого індексу і до молодших індексів
            // використовуємо коли рухаємось від x.Count - 1 до 0

            int j = i - w + 1;
            int pow = 1;
            int temp = 0;
            int max_j = i;
            int max = x[i];

            int count = x.Count - 1;
            while (j < 0)
                j++;

            while (j < i)
            {
                if (x[j] != 0)
                {
                    for (int t = j; t <= i; t++)
                    {
                        temp = temp + pow * x[t];
                        pow = pow << 1;
                    }
                    max_j = j;
                    max = temp;
                    break;
                }
                else
                    j++;
            }
            max_j = i - max_j + 1;

            MyList<int> r = new MyList<int>();
            r.Add(max);
            r.Add(max_j);
            return r;


            /*
            int j = i;
            int pow = 1;
            int temp = x[i];
            int max_j = i;
            int max = x[i];

            while (i - j + 1 < w && j > 0)
            {
                pow = 1;
                temp = 0;
                j--;
                for (int t = j; t <= i; t++)
                {
                    temp = temp + pow * x[t];
                    pow = pow * 2;
                }

                if (temp % 2 != 0)
                {
                    max_j = j;
                    max = temp;
                }
            }
            max_j = i - max_j + 1;

            MyList<int> r = new MyList<int>();
            r.Add(max);
            r.Add(max_j);
            return r;
            */
        }

        public static MyList<BigInteger> NAFRLTable(BigInteger found, BigInteger mod, BigInteger power, int w)
        {
            var table = new MyList<BigInteger>();

            for (int i = 1; i <= power; i += 2)
                table.Add(BinaryRL(found, i, mod));
            return table;
        }
        public static MyList<BigInteger> NAFLRTable(BigInteger found, BigInteger mod, BigInteger power, int w)
        {
            var table = new MyList<BigInteger>();

            for (int i = 1; i <= power; i += 2)
                table.Add(BinaryLR(found, i, mod));
            return table;
        }
        public static MyList<BigInteger> SlidingRLTable(BigInteger found, BigInteger mod, BigInteger power, int w)
        {
            var table = new MyList<BigInteger>();
            table.Add(BinaryRL(found, power, mod));
            for (BigInteger i = 0; i < power - 1; i++)
                table.Add(table[i] * found % mod);
            return table;
        }
        public static MyList<BigInteger> SlidingLRTable(BigInteger found, BigInteger mod, BigInteger power, int w)
        {
            var table = new MyList<BigInteger>();
            table.Add(BinaryLR(found, power, mod));
            for (BigInteger i = 0; i < power - 1; i++)
                table.Add(table[i] * found % mod);
            return table;
        }

        public static List<string> windows(string s, int w)
        {
            string st = s;

            while (st.Length % w != 0)
                st = "0" + st;

            List<string> bins = new List<string>();
            for (int i = 0; i < st.Length / w; i++)
                bins.Add(st.Substring(i * w, w));

            return bins;
        }
        public static MyList<int> ToWNAF(BigInteger power, int w)
        {
            var res = new MyList<int>();
            BigInteger k = power;
            int r;

            int temp;
            int pow2_w = 1 << w;
            int mask = pow2_w - 1;
            int pow2_w_1 = 1 << (w - 1);
            while (k >= 1)
            {
                temp = (int)(k & 1); // k % 2
                if (temp != 0)
                {
                    r = (int)(k & mask); // k % pow2_w

                    if (r >= pow2_w_1)
                        res.Insert(0, r - pow2_w);
                    else
                        res.Insert(0, r);

                    k = k - res[0];
                }
                else
                    res.Insert(0, 0);
                k = k >> 1;
            }
            return res;
        }
        #endregion
    }

    public partial class MainForm : Form
    {
        void Show(BigInteger num, BigInteger pow, BigInteger mod, int window, double table)
        {

            BigInteger nsd = NSD(mod, num % mod);

            bool nsdMethods = false;
            string methodNumber = "";

            for (int i = 0; i < OperationsListTest.CheckedIndices.Count; i++)
                {
                    #region Binary
                    if (OperationsListTest.CheckedIndices[i] == 0) { OperationsResult.Items.Add("Binary RL\t\t: " + (PowFunctions.BinaryRL(num, pow, mod)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 1) { OperationsResult.Items.Add("Binary LR\t\t: " + (PowFunctions.BinaryLR(num, pow, mod)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 10) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("NAF Binary RL\t\t: " + (PowFunctions.NAFBinaryRL(num, pow, mod).ToString())); OperationsResult.Update();}
                    if (OperationsListTest.CheckedIndices[i] == 11) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("NAF Binary LR\t\t: " + (PowFunctions.NAFBinaryLR(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 16) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-1); continue; }
                        OperationsResult.Items.Add("Add Sub RL\t\t: " + (PowFunctions.AddSubRL(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 17) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-1); continue; }
                        OperationsResult.Items.Add("Add Sub LR\t\t: " + (PowFunctions.AddSubLR(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 18) { OperationsResult.Items.Add("Joye_double_and_add\t: " + (PowFunctions.Joye_double_and_add(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 19) { OperationsResult.Items.Add("MontgomeryLadder\t\t: " + (PowFunctions.MontgomeryLadder(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 20) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-1)+".1"; continue; }
                        OperationsResult.Items.Add("DBNS1RL 1\t\t: " + (PowFunctions.DBNS1RL(num, pow, mod, true, AdditionalParameters.A, AdditionalParameters.B).ToString())); OperationsResult.Update();}
                    if (OperationsListTest.CheckedIndices[i] == 21) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-2)+".2"; continue; }
                        OperationsResult.Items.Add("DBNS1RL 2\t\t: " + (PowFunctions.DBNS1RL(num, pow, mod, false, AdditionalParameters.A, AdditionalParameters.B).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 22) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-2)+".1"; continue; }
                        OperationsResult.Items.Add("DBNS1LR 1\t\t: " + (PowFunctions.DBNS1LR(num, pow, mod, true, AdditionalParameters.A, AdditionalParameters.B).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 23) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3)+".2"; continue; }
                        OperationsResult.Items.Add("DBNS1LR 2\t\t: " + (PowFunctions.DBNS1LR(num, pow, mod, false, AdditionalParameters.A, AdditionalParameters.B).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 24) { OperationsResult.Items.Add("DBNS2RL\t\t: " + (PowFunctions.BinaryRL(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 25) { OperationsResult.Items.Add("DBNS2LR\t\t: " + (PowFunctions.BinaryLR(num, pow, mod).ToString())); OperationsResult.Update(); }
                    #endregion

                    #region Window
                    if (OperationsListTest.CheckedIndices[i] == 2) { OperationsResult.Items.Add("Window RL\t\t: " + (PowFunctions.WindowRL(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 3) { OperationsResult.Items.Add("Window RL Dic\t\t: " + (PowFunctions.WindowRL_Dic(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 4) { OperationsResult.Items.Add("Window LR\t\t: " + (PowFunctions.WindowLR(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 5) { OperationsResult.Items.Add("Window LR Dic\t\t: " + (PowFunctions.WindowLR_Dic(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 6) { OperationsResult.Items.Add("Sliding RL\t\t: " + (PowFunctions.SlidingRL(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 7) { OperationsResult.Items.Add("Sliding RL Dic\t\t: " + (PowFunctions.SlidingRL_Dic(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 8) { OperationsResult.Items.Add("Sliding LR\t\t: " + (PowFunctions.SlidingLR(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 9) { OperationsResult.Items.Add("Sliding LR Dic\t\t: " + (PowFunctions.SlidingLR_Dic(num, pow, mod, window, out table)).ToString()); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 12) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("NAF Sliding RL\t\t: " + (PowFunctions.NAFSlidingRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 13) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("NAF Sliding LR\t\t: " + (PowFunctions.NAFSlidingLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 14) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("wNAF RL\t\t: " + (PowFunctions.NAFWindowRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 15) { if (nsd != 1) { nsdMethods = true; methodNumber +=" "+ (OperationsListTest.CheckedIndices[i]-3); continue; }
                        OperationsResult.Items.Add("wNAF LR\t\t: " + (PowFunctions.NAFWindowLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    //if (OperationsListTest.CheckedIndices[i] == 16 && nsd == 1) { if (nsd == 1) { nsdMethods = true; continue; } OperationsResult.Items.Add("wNAF Slide RL\t\t: " + (PowFunctions.wNAFSlidingRL(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    //if (OperationsListTest.CheckedIndices[i] == 17 && nsd == 1) { if (nsd == 1) { nsdMethods = true; continue; } OperationsResult.Items.Add("wNAF Slide LR\t\t: " + (PowFunctions.wNAFSlidingLR(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 26) { OperationsResult.Items.Add("ModWindow LR1\t\t: " + (PowFunctions.WindowLRMod1(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 27) { OperationsResult.Items.Add("ModWindow LR2\t\t: " + (PowFunctions.WindowLRMod2(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 28) { OperationsResult.Items.Add("ModWindow LR3\t\t: " + (PowFunctions.WindowLRMod3(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 29) { OperationsResult.Items.Add("ModWindow LR\t\t: " + (PowFunctions.WindowLRMod(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 30) { OperationsResult.Items.Add("ModWindow LR1(shift)\t: " + (PowFunctions.WindowLRMod1_Shift(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 31) { OperationsResult.Items.Add("ModWindow LR2(shift)\t: " + (PowFunctions.WindowLRMod2_Shift(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 32) { OperationsResult.Items.Add("ModWindow LR3(shift)\t: " + (PowFunctions.WindowLRMod3_Shift(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 33) { OperationsResult.Items.Add("ModWindow LR (shift)\t: " + (PowFunctions.WindowLRMod_Shift(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 34) { OperationsResult.Items.Add("ModWindow LR1(upgrade)\t: " + (PowFunctions.WindowLRMod1_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 35) { OperationsResult.Items.Add("ModWindow LR2(upgrade)\t: " + (PowFunctions.WindowLRMod2_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 36) { OperationsResult.Items.Add("ModWindow LR3(upgrade)\t: " + (PowFunctions.WindowLRMod3_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 37) { OperationsResult.Items.Add("ModWindow LR (upgrade)\t: " + (PowFunctions.WindowLRMod_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 38) { OperationsResult.Items.Add("ModWindow LR1(NoBinary)\t: " + (PowFunctions.WindowLRMod1_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 39) { OperationsResult.Items.Add("ModWindow LR2(NoBinary)\t: " + (PowFunctions.WindowLRMod2_NoBinary(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 40) { OperationsResult.Items.Add("ModWindow LR3(NoBinary)\t: " + (PowFunctions.WindowLRMod3_NoBinary(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 41) { OperationsResult.Items.Add("ModWindow LR (NoBinary)\t: " + (PowFunctions.WindowLRMod_NoBinary(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 42) { OperationsResult.Items.Add("ModWindow LR1(Final)\t: " + (PowFunctions.WindowLRMod1_Upgrade(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 43) { OperationsResult.Items.Add("ModWindow LR2(Final)\t: " + (PowFunctions.WindowLRMod2_Final(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 44) { OperationsResult.Items.Add("ModWindow LR3(Final)\t: " + (PowFunctions.WindowLRMod3_Final(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 45) { OperationsResult.Items.Add("ModWindow LR (Final)\t: " + (PowFunctions.WindowLRMod_Final(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 46) { OperationsResult.Items.Add("Sliding Prime with_11_11\t: " + (PowFunctions.Sliding_Prime(num, pow, mod, window, out table, true).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 47) { OperationsResult.Items.Add("Sliding Prime\t\t: " + (PowFunctions.Sliding_Prime(num, pow, mod, window, out table, false).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 48) { OperationsResult.Items.Add("Adaptive\t\t\t: " + (PowFunctions.Adaptive(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 49) { OperationsResult.Items.Add("Pow C#\t\t\t: " + (BigInteger.ModPow(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 50) { OperationsResult.Items.Add("Bonus1\t\t\t: " + (PowFunctions.Bonus1(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 51) { OperationsResult.Items.Add("Bonus2\t\t\t: " + (PowFunctions.Bonus2(num, pow, mod).ToString())); OperationsResult.Update(); }
                    if (OperationsListTest.CheckedIndices[i] == 52) { OperationsResult.Items.Add("Bonus Window\t\t: " + (PowFunctions.Bonus(num, pow, mod, window, out table).ToString())); OperationsResult.Update(); }

                    #endregion
                }

                if (nsd != 1 && nsdMethods)
                {
                    string message = "Введені значення основи та модуля не є взаємнопростими, а для методу(-ів) "
                                     + methodNumber + " вони мають бути взаємнопростими";
                    MessageBox.Show(message, "Warning", MessageBoxButtons.OK, MessageBoxIcon.Information);
                }
            }
    }
}


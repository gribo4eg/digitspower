using System;
using System.Collections.Generic;
using System.Numerics;
using System.Security.Cryptography;
using System.Diagnostics;
using System.Windows.Forms;
using static DigitsPower.HelpMethods;
using static System.Math;
using static DigitsPower.MontgomeryMethods;

namespace DigitsPower
{
    public delegate BigInteger Inverse(BigInteger mod, BigInteger found);
    public delegate BigInteger Multiply(BigInteger a, BigInteger b, BigInteger m, BigInteger inv);
    public delegate BigInteger OutRes(BigInteger res, BigInteger mod, BigInteger inv);
    public delegate BigInteger InRes(ref BigInteger a, ref BigInteger b, BigInteger m);

    public class MyList<T> : List<T>
    {
        public T this[BigInteger index]    // Indexer declaration
        {
            get { return base[(int)index]; }
            set { base[(int)index] = value; }
        }
    }
    static class PowFunctions
    {
        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

        #region binary

        public static BigInteger mul(BigInteger a, BigInteger b, BigInteger m)
        {
            return ((a * b) % m);
        }

        private static bool GetMontgomery{get { return (Application.OpenForms[0] as MainForm).MontgomeryChecked; } }

        public static BigInteger BinaryRL(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;

            found = found % mod;
            BigInteger t = found;
            BigInteger res = 1;

            string binary = ConvToBinary(pow); // старший біт в рядку під номером 0

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);
                for (int i = binary.Length - 1; i >= 0; i--)
                {
                    if (binary[i] == '1')
                        res = MontgomeryMethods.MontgomeryMultDomain(t, res, mod, parameters);
                    t = MontgomeryMethods.MontgomeryMultDomain(t, t, mod, parameters); ;
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                for (int i = binary.Length - 1; i >= 0; i--)
                {
                    if (binary[i] == '1')
                        res = t * res % mod; //Приведення за модулем після кожного кроку
                    t = t * t % mod;
                }
            }

            return res;

        }

        public static BigInteger BinaryLR(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;
            BigInteger res, t;
            string binary = ConvToBinary(pow); // старший біт в рядку під номером 0

            res = 1;
            t = found % mod;

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);

                for (int i = 0; i < binary.Length; i++)
                {
                    res = MontgomeryMethods.MontgomeryMultDomain(res, res, mod, parameters);
                    if ('1' == binary[i])
                        res = MontgomeryMethods.MontgomeryMultDomain(t, res, mod, parameters);
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                for (int i = 0; i < binary.Length; i++)
                {
                    res = res * res % mod;
                    if (binary[i] == '1')
                        res = t * res % mod;//Приведення до степеня після кожного кроку
                }
            }

            return res;
        }


        public static BigInteger NAFBinaryRL(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;

            BigInteger res, t;
            MyList<int> x;
            found = found % mod;

            res = 1;
            t = found;
            x = ToNAF(pow);

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);

                for (int i = 0; i < x.Count; i++)
                {
                    if (x[i] == 1)
                        res = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                    else if (x[i] == -1)
                        res = MontgomeryMethods.MontgomeryMultDomain(res, MontgomeryMethods.MontgomeryInverse(mod, t, parameters), mod, parameters);
                    t = MontgomeryMethods.MontgomeryMultDomain(t, t, mod, parameters);
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                for (int i = 0; i < x.Count; i++)
                {
                    if (x[i] == 1)
                        res = res * t % mod;
                    else if (x[i] == -1)
                        res = res * Euclid_2_1(mod, t) % mod;
                    t = t * t % mod;
                }
            }


            return res;

        }

        public static BigInteger NAFBinaryLR(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;

            BigInteger res;
            MyList<int> x;
            found = found % mod;

            x = ToNAF(pow);
            if (mont)
            {
                BigInteger t;

                res = 1;
                t = found;
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);
                for (int i = x.Count - 1; i >= 0; i--)
                {
                    res = MontgomeryMethods.MontgomeryMultDomain(res, res, mod, parameters);
                    if (x[i] == 1)
                        res = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                    else if (x[i] == -1)
                        res = MontgomeryMethods.MontgomeryMultDomain(res, MontgomeryInverse(mod, t, parameters), mod, parameters);
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                BigInteger inv = Euclid_2_1(mod, found);
                res = found;

                for (int i = x.Count - 2; i > -1; i--)
                {
                    res = res * res % mod;
                    if (x[i] == 1)
                        res = res * found % mod;
                    else if (x[i] == -1)
                        res = res * inv % mod;
                }
            }
            return res;

        }

        public static BigInteger AddSubRL(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;
            found = found % mod;

            BigInteger res, pow3;
            string pow3_bin, pow_bin;

            res = 1;
            pow3 = pow * 3;

            pow_bin = ConvToBinary(pow);
            pow3_bin = ConvToBinary(pow3);
            int pow3_length = pow3_bin.Length - 1;
            while (pow_bin.Length <= pow3_length)
                pow_bin = '0' + pow_bin;

            BigInteger t = found;
            if (mont)
            {
                MyList<BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);

                for(int i = pow3_length - 1; i >= 0; i--)
                {
                    if (pow3_bin[i] == '1' && pow_bin[i] == '0')
                        res = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                    else if (pow3_bin[i] == '0' && pow_bin[i] == '1')
                        res = MontgomeryMethods.MontgomeryMultDomain(res, MontgomeryInverse(mod, t, parameters), mod, parameters);
                    t = MontgomeryMethods.MontgomeryMultDomain(t, t, mod, parameters);
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                for (int i = pow3_length - 1; i > -1; i--)
                {
                    if (pow3_bin[i] == '1' && pow_bin[i] == '0')
                        res = res * t % mod;
                    else if (pow3_bin[i] == '0' && pow_bin[i] == '1')
                        res = res * Euclid_2_1(mod, t) % mod;
                    t = t * t % mod;
                }
            }
            return res;
        }

        public static BigInteger AddSubLR(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;
            found = found % mod;

            BigInteger res, pow3;
            string pow3_bin, pow_bin;

            res = found;
            pow3 = pow * 3;

            pow_bin = ConvToBinary(pow);
            pow3_bin = ConvToBinary(pow3);
            int pow3_length = pow3_bin.Length - 1;
            while (pow_bin.Length <= pow3_length)
                pow_bin = '0' + pow_bin;

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref found, ref res, mod);

                for (int i = 1; i < pow3_length; i++)
                {
                    res = MontgomeryMultDomain(res, res, mod, parameters);
                    if (pow3_bin[i] == '1' && pow_bin[i] == '0')
                        res = MontgomeryMultDomain(res, found, mod, parameters);
                    else if (pow3_bin[i] == '0' && pow_bin[i] == '1')
                        res = MontgomeryMultDomain(res, MontgomeryInverse(mod, found, parameters), mod, parameters);
                }

                res = outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                BigInteger inv = Euclid_2_1(mod, found);

                for (int i = 1; i < pow3_length; i++)
                {
                    res = res * res % mod;
                    if (pow3_bin[i] == '1' && pow_bin[i] == '0')
                        res = res * found % mod;
                    else if (pow3_bin[i] == '0' && pow_bin[i] == '1')
                        res = res * inv % mod;
                }
            }

            return res;
        }

        public static BigInteger Joye_double_and_add(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;

            BigInteger res, t;
            int powLen;
            string pow_bin;
            found = found % mod;

            pow_bin = ConvToBinary(pow);
            res = 1;
            t = found;
            if (mont)
            {

                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);

                powLen = (int) (Log((double) pow, 2) + 1);
                for (int i = 0; i < powLen; i++)
                {
                    if ('1' == pow_bin[powLen - (int) i - 1])
                    {
                        res = MontgomeryMethods.MontgomeryMultDomain(res, res, mod, parameters);
                        res = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                    }
                    else
                    {
                        t = MontgomeryMethods.MontgomeryMultDomain(t, t, mod, parameters);
                        t = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                    }
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                powLen = pow_bin.Length;
                for (int i = powLen - 1; i >= 0; i--)
                {
                    if (pow_bin[i] == '1')
                    {
                        res = res * res % mod;
                        res = res * t % mod;
                    }
                    else
                    {
                        t = t * t % mod;
                        t = res * t % mod;
                    }
                }
            }

            return res;
        }

        public static BigInteger MontgomeryLadder(BigInteger found, BigInteger pow, BigInteger mod)
        {
            bool mont = GetMontgomery;
            BigInteger res, t;
            int powLen;
            string pow_bin;
            found = found % mod;

            pow_bin = ConvToBinary(pow);
            t = found;
            res = 1;
            if(mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref t, ref res, mod);

                powLen = (int)(Log((double)pow, 2));
                for (int i = powLen; 0 <= i; i--)
                {
                    if (pow_bin[powLen - (int)i] == '0')
                    {
                        t = MontgomeryMethods.MontgomeryMultDomain(t, res, mod, parameters);
                        res = MontgomeryMethods.MontgomeryMultDomain(res, res, mod, parameters);
                    }
                    else
                    {
                        res = MontgomeryMethods.MontgomeryMultDomain(res, t, mod, parameters);
                        t = MontgomeryMethods.MontgomeryMultDomain(t, t, mod, parameters);
                    }
                }
                res = MontgomeryMethods.outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                powLen = pow_bin.Length;
                for (int i = 0; i < powLen; i++)
                {
                    if (pow_bin[i] == '0')
                    {
                        t = t * res % mod;
                        res = res * res % mod;
                    }
                    else
                    {
                        res = res * t % mod;
                        t = t * t % mod;
                    }
                }
            }

            return res;
        }

        public static BigInteger DBNS1RL(BigInteger found, BigInteger pow, BigInteger mod, bool convert_method, long A, long B)
        {
            BigInteger[,] mas_k;
            mas_k = convert_method ? Convert_to_DBNS_1(pow, A, B)
                                   : Convert_to_DBNS_2(pow, A, B);
            found = found % mod;

            long lastindex = mas_k.GetLength(0) - 1;
            BigInteger t = found;
            BigInteger res = 1;

            for (long i = 0; i < mas_k[lastindex, 1]; i++)
                t = mul(t, t, mod);


            for (long i = 0; i < mas_k[lastindex, 2]; i++)
                t = mul(mul(t, t, mod), t, mod);

            if (mas_k[lastindex, 0] == -1)
                res = mul(res, Euclid_2_1(mod, t), mod);
            else if (mas_k[lastindex, 0] == 1)
                res = mul(res, t, mod);

            for (long i = lastindex - 1; i >= 0; i--)
            {
                BigInteger u = mas_k[i, 1] - mas_k[i + 1, 1];
                BigInteger v = mas_k[i, 2] - mas_k[i + 1, 2];
                for (long j = 0; j < u; j++)
                    t = mul(t, t, mod);

                for (long j = 0; j < v; j++)
                    t = mul(mul(t, t, mod), t, mod);

                if (mas_k[i, 0] == -1)
                    res = mul(res, Euclid_2_1(mod, t), mod);
                else if (mas_k[i, 0] == 1)
                    res = mul(res, t, mod);
            }
            return res;
        }

        public static BigInteger DBNS1LR(BigInteger found, BigInteger pow, BigInteger mod, bool convert_method, long A, long B)
        {
            BigInteger[,] mas_k;
            found = found % mod;

            BigInteger t = found;
            BigInteger res = 1;

            mas_k = convert_method ? Convert_to_DBNS_1(pow, A, B)
                                   : Convert_to_DBNS_2(pow, A, B);

            if (mas_k[0, 0] == -1)
                res = Euclid_2_1(mod, t);
            else if (mas_k[0, 0] == 1)
                res = t;

            for (long i = 0; i < mas_k.GetLength(0) - 1; i++)
            {
                BigInteger u = mas_k[i, 1] - mas_k[i + 1, 1];
                BigInteger v = mas_k[i, 2] - mas_k[i + 1, 2];
                for (long j = 0; j < u; j++)
                    res = mul(res, res, mod);

                for (long j = 0; j < v; j++)
                    res = mul(mul(res, res, mod), res, mod);

                if (mas_k[i + 1, 0] < 0)
                    res = mul(res, Euclid_2_1(mod, t), mod);
                else
                    res = mul(res, t, mod);
            }

            for (long i = 0; i < mas_k[mas_k.GetLength(0) - 1, 1]; i++)
                res = mul(res, res, mod);

            for (long i = 0; i < mas_k[mas_k.GetLength(0) - 1, 2]; i++)
                res = mul(mul(res, res, mod), res, mod);
            return res;
        }

        public static BigInteger DBNS2RL(BigInteger found, BigInteger pow, BigInteger mod)
        {
            MyList<int[]> mas_k = ToDBNS2RL(pow);
            found = found % mod;

            BigInteger res = 1;
            BigInteger t = found;

            for (int i = 0; i < mas_k.Count; i++)
            {
                for (int j = 0; j < mas_k[i][1]; j++)
                    t = mul(t, t, mod);
                for (int j = 0; j < mas_k[i][2]; j++)
                    t = mul((mul(t, t, mod)), t, mod);

                if (mas_k[i][0] == 1)
                    res = mul(res, t, mod);
                else if (mas_k[i][0] == -1)
                    res = mul(res, Euclid_2_1(t, mod), mod);
            }
            return res;
        }

        public static BigInteger DBNS2LR(BigInteger found, BigInteger pow, BigInteger mod)
        {
            MyList<int[]> mas_k = ToDBNS2LR(pow);
            found = found % mod;
            BigInteger res = found;

            for (int i = mas_k.Count - 1; i > 0; i--)
            {
                for (int j = 0; j < mas_k[i][1]; j++)
                    res = mul(res, res, mod);
                for (int j = 0; j < mas_k[i][2]; j++)
                    res = mul(mul(res,  res, mod), res, mod);

                if (mas_k[i][0] == 1)
                    res = mul(res, found, mod);
                else if (mas_k[i][0] == -1)
                    res = mul(res, Euclid_2_1(found, mod), mod);
            }
            for (int j = 0; j < mas_k[0][1]; j++)
                res = mul(res, res, mod);
            for (int j = 0; j < mas_k[0][2]; j++)
                res = mul(mul(res,  res, mod), res, mod);
            return res;
        }

        public static BigInteger Bonus1(BigInteger found, BigInteger pow, BigInteger mod)
        {
            BigInteger res, t, inverse;
            int powLen;
            string pow_bin;
            found = found % mod;

            pow_bin = ConvToBinary(pow);
            res = 1;
            t = found;
            inverse = AdditionalParameters.inRes(ref t, ref res, mod);
            powLen = (int)(Log((double)pow, 2));
            for (int i = powLen; 0 <= i; i--)
            {
                res = AdditionalParameters.mul(res, res, mod, inverse);
                if (pow_bin[powLen - (int)i] == '1')
                    res = AdditionalParameters.mul(t, res, mod, inverse);
            }
            res = AdditionalParameters.outRes(res, mod, inverse);
            return res;
        }
        public static BigInteger Bonus2(BigInteger found, BigInteger pow, BigInteger mod)
        {
            BigInteger res, t, inverse;
            int powLen;
            string pow_bin;
            found = found % mod;

            pow_bin = ConvToBinary(pow);

            res = 1;
            t = found;
            inverse = AdditionalParameters.inRes(ref t, ref res, mod);
            powLen = (int)(Log((double)pow, 2));
            for (int i = powLen; 0 <= i; i--)
            {
                res = AdditionalParameters.mul(res, res, mod, inverse);
                if (pow_bin[powLen - (int)i] == '1')
                    res = AdditionalParameters.mul(t, res, mod, inverse);
            }
            res = AdditionalParameters.outRes(res, mod, inverse);
            return res;
        }
        #endregion
        #region window
        private static MyList<BigInteger> Table(BigInteger found, BigInteger pow, int w, BigInteger mod)
        {
            var table = new MyList<BigInteger>();

            table.Add(found);
            for (BigInteger i = 0; i < BigInteger.Parse((Pow(2, w) - 2).ToString()); i++)
                table.Add((table[i] * found) % mod);
            return table;
        }

        public static BigInteger WindowRL(BigInteger found, BigInteger pow,  BigInteger mod, int w, out double table_time)
        {
            bool mont = GetMontgomery;
            // 3.1
            found = found % mod;

            int i;
            int count_elem = 1 << w;
            MyList<BigInteger> table = new MyList<BigInteger>();

            BigInteger a = found;
            BigInteger res = 1;

            List<string> bins = windows(ConvToBinary(pow), w);

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref a, ref res, mod);

                for (i = 1; i <= count_elem; i++)
                    table.Add(res);

                int index;
                for (i = bins.Count - 1; i > -1; i--)
                {
                    index = Convert.ToInt32(bins[i], 2);

                    if (index != 0)
                        table[index - 1] = MontgomeryMultDomain(table[index - 1], a, mod, parameters);

                    for (int k = 0; k < w; k++)
                        a = MontgomeryMultDomain(a, a, mod, parameters);
                }

                Stopwatch stw = new Stopwatch();
                stw.Start();
                for (i = count_elem - 1; i > 0; i--)
                {
                    table[i - 1] = MontgomeryMultDomain(table[i - 1], table[i], mod, parameters);
                    table[0] = MontgomeryMultDomain(table[0], table[i], mod, parameters);
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                res = outMontgomeryDomain(table[0], mod, parameters);
                return res;
            }
            else
            {
                for (i = 1; i <= count_elem; i++)
                    table.Add(1);

                int index;
                for (i = bins.Count - 1; i > -1; i--)
                {
                    index = Convert.ToInt32(bins[i], 2);

                    if (index != 0)
                        table[index - 1] = table[index - 1] * a % mod;

                    for (int k = 0; k < w; k++)
                        a = a * a % mod;
                }

                Stopwatch stw = new Stopwatch();
                stw.Start();
                for (i = count_elem - 1; i > 0; i--)
                {
                    table[i - 1] = table[i - 1] * table[i] % mod;
                    table[0] = table[0] * table[i] % mod;
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                return table[0];
            }

            /*
            // залишити закоментованим
            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            MyList<BigInteger> table = Table(found, pow, w, mod);
            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;

            BigInteger res = 1;

            List<string> bins = windows(ConvToBinary(pow), w);

            for (int i = bins.Count - 1; i > -1; i--)
            {
                int c = Convert.ToInt32(bins[i], 2);
                if (c != 0)
                    res = mul(res,  table[c - 1], mod);

                for (int k = 0; k < w; k++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = mul(table[j], table[j], mod);
            }
            return res;
            */
        }

        public static BigInteger WindowRL_Dic(BigInteger found, BigInteger pow,  BigInteger mod, int w, out double table_time)
        {
            // 3.2. Dic
            bool mont = GetMontgomery;
            found = found % mod;

            int i;
            int count_elem = 1 << w;
            Dictionary<string, BigInteger> table = new Dictionary<string, BigInteger>();

            string temp;
            BigInteger res = 1;
            BigInteger a = found;

            string pow_bin = ConvToBinary(pow);
            string str_null = "";
            string temp_previous;
            string str_one = "";

            Stopwatch stw = new Stopwatch();

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref a, ref res, mod);

                for (i = 1; i < count_elem; i++)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    table.Add(temp, res);
                }

                while (pow_bin.Length % w != 0)
                    pow_bin = "0" + pow_bin;

                for (i = 0; i < w; i++)
                    str_null = "0" + str_null;


                for (i = pow_bin.Length - w; i >= 0; i = i - w)
                {
                    string subs = pow_bin.Substring(i, w);
                    if (subs != str_null) table[subs] = MontgomeryMultDomain(table[subs], a, mod, parameters);

                    for (int k = 0; k < w; k++)
                        a = MontgomeryMultDomain(a, a, mod, parameters);
                }

                stw.Start();

                for (i = 1; i < w; i++)
                    str_one = "0" + str_one;
                str_one = str_one + "1";

                for (i = count_elem - 1; i > 1; i--)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    temp_previous = ConvToBinary(i - 1);
                    while (temp_previous.Length != w)
                        temp_previous = "0" + temp_previous;
                    table[temp_previous] = MontgomeryMultDomain(table[temp_previous], table[temp], mod, parameters);
                    table[str_one] = MontgomeryMultDomain(table[str_one], table[temp], mod, parameters);
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                res = outMontgomeryDomain(table[str_one], mod, parameters);
                return res;
            }
            else
            {
                for (i = 1; i < count_elem; i++)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    table.Add(temp, 1);
                }

                while (pow_bin.Length % w != 0)
                    pow_bin = "0" + pow_bin;

                for (i = 0; i < w; i++)
                    str_null = "0" + str_null;


                for (i = pow_bin.Length - w; i >= 0; i = i - w)
                {
                    string subs = pow_bin.Substring(i, w);
                    if (subs != str_null) table[subs] = table[subs] * a % mod;

                    for (int k = 0; k < w; k++)
                        a = a * a % mod;
                }

                stw.Start();

                for (i = 1; i < w; i++)
                    str_one = "0" + str_one;
                str_one = str_one + "1";

                for (i = count_elem - 1; i > 1; i--)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    temp_previous = ConvToBinary(i - 1);
                    while (temp_previous.Length != w)
                        temp_previous = "0" + temp_previous;
                    table[temp_previous] = table[temp_previous] * table[temp] % mod;
                    table[str_one] = table[str_one] * table[temp] % mod;
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                return table[str_one];
            }
        }

        public static BigInteger WindowLR(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            bool mont = GetMontgomery;
            Stopwatch stw = new Stopwatch();
            BigInteger res = 1;

            found = found % mod;
            MyList<BigInteger> table = new MyList<BigInteger>();
            List<string> bins = windows(ConvToBinary(pow), w);

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref found, ref res, mod);

                stw.Start();

                table.Add(found);

                int count = (1 << w) - 2;
                for (BigInteger i = 0; i < count; i++)
                    table.Add(MontgomeryMultDomain(table[i], found, mod, parameters));

                stw.Stop();

                table_time = stw.Elapsed.TotalMilliseconds;

                for (int i = 0; i < bins.Count; i++)
                {
                    for (int k = 0; k < w; k++)
                        res = MontgomeryMultDomain(res, res, mod, parameters);

                    int c = Convert.ToInt32(bins[i], 2);
                    if (c != 0) res = MontgomeryMultDomain(res, table[c - 1], mod, parameters);
                }

                res = outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                stw.Start();

                table.Add(found);

                int count = (1 << w) - 2;
                for (BigInteger i = 0; i < count; i++)
                    table.Add((table[i] * found) % mod);

                stw.Stop();

                table_time = stw.Elapsed.TotalMilliseconds;

                for (int i = 0; i < bins.Count; i++)
                {
                    for (int k = 0; k < w; k++)
                        res = mul(res, res, mod);

                    int c = Convert.ToInt32(bins[i], 2);
                    if (c != 0) res = mul(res, table[c - 1], mod);
                }
            }
            return res;
        }

        public static BigInteger WindowLR_Dic(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            bool mont = GetMontgomery;
            Stopwatch stw = new Stopwatch();
            found = found % mod;

            BigInteger res = 1;

            string temp;
            BigInteger temp_value;
            int temp_size = (1 << w) - 1;

            string pow_bin = ConvToBinary(pow);
            string str_null = "";
            Dictionary < string, BigInteger> table = new Dictionary<string, BigInteger>();
            temp = ConvToBinary(1);

            if(mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref found, ref res, mod);

                stw.Start();

                while (temp.Length != w)
                    temp = "0" + temp;

                table.Add(temp, found);
                temp_value = found;

                for (int i = 2; i <= temp_size; i++)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    temp_value = MontgomeryMultDomain(temp_value, found, mod, parameters);
                    table.Add(temp, temp_value);
                }

                stw.Stop();

                table_time = stw.Elapsed.TotalMilliseconds;

                while (pow_bin.Length % w != 0)
                    pow_bin = "0" + pow_bin;

                for (int i = 0; i < w; i++)
                    str_null = "0" + str_null;

                for (int i = 0; i < pow_bin.Length; i = i + w)
                {
                    for (int k = 0; k < w; k++)
                        res = MontgomeryMultDomain(res, res, mod, parameters);

                    string subs = pow_bin.Substring(i, w);
                    if (subs != str_null) res = MontgomeryMultDomain(res, table[subs], mod, parameters);
                }

                res = outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                stw.Start();

                while (temp.Length != w)
                    temp = "0" + temp;

                table.Add(temp, found);
                temp_value = found;

                for (int i = 2; i <= temp_size; i++)
                {
                    temp = ConvToBinary(i);
                    while (temp.Length != w)
                        temp = "0" + temp;
                    temp_value = temp_value * found % mod;
                    table.Add(temp, temp_value);
                }

                stw.Stop();

                table_time = stw.Elapsed.TotalMilliseconds;

                while (pow_bin.Length % w != 0)
                    pow_bin = "0" + pow_bin;

                for (int i = 0; i < w; i++)
                    str_null = "0" + str_null;

                for (int i = 0; i < pow_bin.Length; i = i + w)
                {
                    for (int k = 0; k < w; k++)
                        res = res * res % mod;

                    string subs = pow_bin.Substring(i, w);
                    if (subs != str_null) res = res * table[subs] % mod;
                }
            }
            return res;

        }

        public static BigInteger WindowLRMod1(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);
            found = found % mod;

            stw = new Stopwatch();
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 1; i < w; i++)
                table.Add(mul(table[i - 1], table[i - 1], mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                res = mul(res, res, mod);
                if ('1' == pow_bin[powLen - 1 - (int)i])
                {
                    if (0 < i && '0' == pow_bin[powLen - (int)i - 1])
                    {
                        t = 1;
                        while (t < w
                            && t < i
                            && '0' == pow_bin[powLen - (int)i - t])
                        {
                            t++;
                        }
                        res = BinaryRL(res, 1 << (t - 1), mod);
                        res = mul(res, table[t - 1], mod);
                        i = i - t + 1;
                    }
                    else
                    {
                        res = mul(res, found, mod);
                    }
                }
                --i;
            }
            return res;
        }

        public static BigInteger WindowLRMod2(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 1; i <= w; i++)
                table.Add(BinaryRL(found, (1 << i) - 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                res = mul(res, res, mod);
                if ('1' == pow_bin[powLen - 1 - (int)i])
                {
                    t = 1;
                    while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                        t++;
                    res = BinaryRL(res, 1 << (t - 1), mod);
                    res = mul(res, table[t], mod);
                    i = i - t + 1;
                }

                    i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod3(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 2; i <= w; i++)
                table.Add(BinaryRL(found, (1 << (i - 1)) + 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)//>=
            {
                res = mul(res, res, mod);
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            res = BinaryRL(res, 1 << (t - 1), mod);
                            res = mul(res, table[t - 1], mod);
                            i = i - t + 1;
                        }
                        else
                            res = mul(res, table[0], mod);

                    }
                    else
                    {
                        res = mul(res, res, mod);
                        res = mul(res, table[1], mod);
                        i--;
                    }

                }
                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            string pow_bin;
            BigInteger[] table1, table2;

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            table1 = new BigInteger[w+1];
            table2 = new BigInteger[w+1];
            table1[0] = table2[0] = found;
            for (i = 1; i < w; i++) {
                table1[i] = BinaryRL(found, (1 << (i + 1)) - 1, mod);
                table2[i] = BinaryRL(found, (1 << i) + 1, mod);
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            pow_bin = ConvToBinary(pow);
            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            res = BinaryRL(res, 1 << t, mod);
                            res = mul(res, table2[t - 1], mod);
                            i = i - t;
                        }
                        else
                        {
                            res = mul(res, res, mod);
                            res = mul(res, table2[0], mod);
                            i--;
                        }
                    }
                    else
                    {
                        t = 1;
                        while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                            t++;
                        res = BinaryRL(res, 1 << t, mod);
                        res = mul(res, table1[t - 1], mod);
                        i = i - t;
                    }
                }
                else
                {
                    res = mul(res, res, mod);
                    i--;
                }

            }
            return res;
        }

        public static BigInteger WindowLRMod1_Shift(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            var pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 1; i < w; i++)
                table.Add(mul(table[i - 1], table[i - 1], mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Math.Log((double)pow, 2));
            while (i >= 0)
            {
                res = mul(res, res, mod);
                if (1 == ((pow >> i) & 1))
                {
                    if (0 < i && 0 == ((pow >> (i - 1)) & 1))
                    {
                        t = 1;
                        while (t < w
                            && t < i
                            && 0 == ((pow >> (i - t)) & 1))
                        {
                            t++;
                        }
                        res = BinaryRL(res, 1 << (t - 1), mod);
                        res = mul(res, table[t - 1], mod);
                        i = i - t + 1;
                    }
                    else
                    {
                        res = mul(res, found, mod);
                        //--i;
                    }
                }
                --i;
            }
            return res;
        }

        public static BigInteger WindowLRMod2_Shift(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 1; i <= w; i++)
                table.Add(BinaryRL(found, (1 << i) - 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            while (i >= 0)
            {
                res = mul(res, res, mod);
                if (1 == ((pow >> i) & 1))
                {
                    t = 1;
                    while (t < w && t < i && 1 == ((pow >> (i - t)) & 1))
                        t++;
                    res = BinaryRL(res, 1 << (t - 1), mod);
                    res = mul(res, table[t], mod);
                    i = i - t + 1;
                }

                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod3_Shift(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            var pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();

            table.Add(found);
            for (i = 2; i <= w; i++)
                table.Add(BinaryRL(found, (1 << (i - 1)) + 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            while (i >= 0)//>=
            {
                res = mul(res, res, mod);
                if (1 == ((pow >> i) & 1))
                {
                    if (0 == ((pow >> (i - 1)) & 1))
                    {
                        t = 2;
                        if (t < i)
                        {
                            while (t < w && 0 == ((pow >> (i - t)) & 1))
                                t++;
                            if (t < w)
                                t++;
                            else
                                t = 0;

                        }
                        else
                        {
                            t = 0;
                        }

                        if (t > 0)
                        {
                            res = BinaryRL(res, (1 <<(t - 1)), mod);
                            res = mul(res, table[t - 1], mod);
                            i = i - t + 1;
                        }
                        else
                        {
                            res = mul(res, table[0], mod);
                        }

                    }
                    else
                    {
                        res = mul(res, res, mod);
                        res = mul(res, table[1], mod);
                        i--;
                    }

                }
                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod_Shift(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i;
            int t;
            Stopwatch stw;
            BigInteger res;
            string pow_bin;
            BigInteger[] table1, table2;

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            table1 = new BigInteger[w+1];
            table2 = new BigInteger[w+1];
            table1[0] = table2[0] = found;
            for (i = 1; i < w; i++) {
                table1[i] = BinaryRL(found, (1 << (i + 1)) - 1, mod);
                table2[i] = BinaryRL(found, (1 << i) + 1, mod);
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            pow_bin = ConvToBinary(pow);
            res = 1;
            i = (int)(Log((double)pow, 2));
            while (i >= 0)
            {
                if (1 == ((pow >> i) & 1))
                {
                    if (0 == ((pow >> (i-1)) & 1))
                    {
                        t = 2;
                        if (t < i)
                        {
                            while (t < w && 0 == ((pow >> (i - t)) & 1))
                                t++;
                            if (t < w)
                                t++;
                            else
                                t = 0;

                        }
                        else
                        {
                            t = 0;
                        }

                        if (t > 0)
                        {
                            res = BinaryRL(res, 1 << t, mod);
                            res = mul(res, table2[t - 1], mod);
                            i = i - t;
                        }
                        else
                        {
                            res = mul(res, res, mod);
                            res = mul(res, table2[0], mod);
                            i--;
                        }
                    }
                    else
                    {
                        t = 1;
                        while (t < w && t < i && 1 == ((pow >> (i - t)) & 1))
                            t++;
                        res = BinaryRL(res, 1 << t, mod);
                        res = mul(res, table1[t - 1], mod);
                        i = i - t;
                    }
                }
                else
                {
                    res = mul(res, res, mod);
                    i--;
                }

            }
            return res;
        }

        public static BigInteger WindowLRMod1_Upgrade(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            var pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();
            table.Add(found);
            for (i = 1; i < w; i++)
                table.Add((table[i - 1] * table[i - 1]) % mod);

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Math.Log((double)pow, 2));
	        powLen = i + 1;
            while (i >= 0)
            {
                res = (res * res) % mod;
                if ('1' == pow_bin[powLen - 1 - (int)i])
                {
                    if (0 < i && '0' == pow_bin[powLen - (int)i - 1])
                    {
                        t = 1;
                        while (t < w
                            && t < i
                            && '0' == pow_bin[powLen - (int)i - t])
                        {
                            t++;
                        }
                        for (int j = 1; j <= t; j++)
                            res = res * res % mod;
                        res = (res * table[t - 1]) % mod;
                        i = i - t + 1;
                    }
                    else
                    {
                        res = (res * found) % mod;
                    }
                }
                --i;
            }
            return res;
        }

        public static BigInteger WindowLRMod2_Upgrade(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();
            table.Add(found);
            for (i = 1; i <= w; i++)
                table.Add(BinaryRL(found, (1 << i) - 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                res = (res * res) % mod;
                if ('1' == pow_bin[powLen - 1 - (int)i])
                {
                    t = 1;
                    while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                        t++;
                    for (int j = 1; j <= t - 1; j++)
                        res = res * res % mod;
                    res = (res * table[t]) % mod;
                    i = i - t + 1;
                }

                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod3_Upgrade(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table = new List<BigInteger>();
            table.Add(found);
            for (i = 2; i <= w; i++)
                table.Add(BinaryRL(found, (1 << (i - 1)) + 1, mod));
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)//>=
            {
                res = res * res % mod;
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            for (int j = 1; j <= t - 1; j++)
                                res = res * res % mod;
                            res = res * table[t - 1] % mod;
                            i = i - t + 1;
                        }
                        else
                            res = res * table[0] % mod;

                    }
                    else
                    {
                        res = res * res % mod;
                        res = res * table[1] % mod;
                        i--;
                    }

                }
                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod_Upgrade(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res;
            string pow_bin;
            BigInteger[] table1, table2;

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table1 = new BigInteger[w+1];
            table2 = new BigInteger[w+1];
            table1[0] = table2[0] = found;
            for (i = 1; i < w; i++) {
                table1[i] = BinaryRL(found, (1 << (i + 1)) - 1, mod);
                table2[i] = BinaryRL(found, (1 << i) + 1, mod);
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            pow_bin = ConvToBinary(pow);
            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            for (int j = 1; j <= t; j++)
                                res = res * res % mod;
                            res = res * table2[t - 1] % mod;
                            i = i - t;
                        }
                        else
                        {
                            res = res * res % mod;
                            res = res * table2[0] % mod;
                            i--;
                        }
                    }
                    else
                    {
                        t = 1;
                        while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                            t++;
                        for (int j = 1; j <= t; j++)
                            res = res * res % mod;
                        res = res * table1[t - 1] % mod;
                        i = i - t;
                    }
                }
                else
                {
                    res = res * res % mod;
                    i--;
                }

            }
            return res;
        }

        public static BigInteger WindowLRMod2_NoBinary(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res, temp;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            temp = found;
            stw.Start();
            table = new List<BigInteger>();
            table.Add(found);
            table.Add(found);
            for (i = 2; i <= w; i++)
            {
                found = found * found % mod;
                found = found * temp % mod;
                table.Add(found);
            }


            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                res = (res * res) % mod;
                if ('1' == pow_bin[powLen - 1 - (int)i])
                {
                    t = 1;
                    while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                        t++;
                    for (int j = 1; j <= t - 1; j++)
                        res = res * res % mod;
                    res = (res * table[t]) % mod;
                    i = i - t + 1;
                }

                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod3_NoBinary(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, t, powLen;
            Stopwatch stw;
            BigInteger res, temp;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            table = new List<BigInteger>();
            table.Add(found);
            temp = found;
            for (i = 2; i <= w; i++)
            {
                found = found * found % mod;
                table.Add(found * temp % mod);
            }

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)//>=
            {
                res = res * res % mod;
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            for (int j = 1; j <= t - 1; j++)
                                res = res * res % mod;
                            res = res * table[t - 1] % mod;
                            i = i - t + 1;
                        }
                        else
                            res = res * table[0] % mod;

                    }
                    else
                    {
                        res = res * res % mod;
                        res = res * table[1] % mod;
                        i--;
                    }

                }
                i--;
            }
            return res;
        }

        public static BigInteger WindowLRMod_NoBinary(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, j, t, powLen;
            Stopwatch stw;
            BigInteger res, temp, found2;
            string pow_bin;
            BigInteger[] table1, table2;

            stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            table1 = new BigInteger[w+1];
            table2 = new BigInteger[w+1];
            temp = table1[0] = table2[0] = found;
            found2 = found;
            for (i = 1; i < w; i++)
            {
                found = found * found % mod;
                found = found * temp % mod;

                table1[i] = found;

                found2 = found2 * found2 % mod;

                table2[i] = found2 * temp % mod;

            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            pow_bin = ConvToBinary(pow);
            res = 1;
            i = (int)(Log((double)pow, 2));
            powLen = i + 1;
            while (i >= 0)
            {
                if ('1' == pow_bin[powLen - i - 1])
                {
                    if ((i == 0) || (i > 0 && '0' == pow_bin[powLen - i ]))
                    {
                        t = 2;
                        while (t < w && t <= i && '0' == pow_bin[powLen - i + t - 1])
                            t++;
                        if (t < w && t <= i)
                            t++;
                        else
                            t = 0;

                        if (t > 0)
                        {
                            for (j = 1; j <= t; j++)
                                res = res * res % mod;
                            res = res * table2[t - 1] % mod;
                            i = i - t;
                        }
                        else
                        {
                            res = res * res % mod;
                            res = res * table2[0] % mod;
                            i--;
                        }
                    }
                    else
                    {
                        t = 1;
                        while (t < w && t < i && '1' == pow_bin[powLen - (int)i + t - 1])
                            t++;

                        for (j = 1; j <= t; j++)
                            res = res * res % mod;
                        res = res * table1[t - 1] % mod;
                        i = i - t;
                    }
                }
                else
                {
                    res = res * res % mod;
                    i--;
                }

            }
            return res;
        }

        public static BigInteger WindowLRMod2_Final(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, j, t, powLen;
            Stopwatch stw;
            BigInteger res, temp;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            found = found % mod;
            stw = new Stopwatch();
            temp = found;
            stw.Start();
            table = new List<BigInteger>();
            table.Add(found);
            for (i = 2; i <= w; i++)  // нульовий елемент саме число, далі 11, 111, ...
            {
                found = found * found % mod;
                found = found * temp % mod;
                table.Add(found);
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            powLen = pow_bin.Length;
            int ind;
            i = 0;
            while (i < powLen)
            {
                res = res * res % mod;
                if (pow_bin[i] == '1')
                {
                    t = 1;  // кількість одиниць
                    ind = i + t;
                    while (t < w && ind < powLen && pow_bin[ind] == '1')
                    {
                        t++;
                        ind++;
                    }

                    for (j = 1; j < t; j++)
                        res = res * res % mod;
                    res = res * table[t - 1] % mod;
                    i = i + t;
                }
                else
                    i++;
            }
            return res;
        }

        public static BigInteger WindowLRMod3_Final(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, j, t, powLen;
            Stopwatch stw;
            BigInteger res, temp;
            List<BigInteger> table;
            string pow_bin = ConvToBinary(pow);

            found = found % mod;
            stw = new Stopwatch();
            stw.Start();

            table = new List<BigInteger>();
            table.Add(found);
            temp = found;
            for (i = 2; i <= w; i++)   // нульовий елемент саме число, далі 11, 101, 1001, ...
            {
                found = found * found % mod;
                table.Add(found * temp % mod);
            }

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            res = 1;
            i = 0;
            powLen = pow_bin.Length;
            int ind;
            int count_zero_limit = w - 2;
            while (i < powLen)
            {
                res = res * res % mod;
                if (pow_bin[i] == '1')
                {
                    i++;
                    t = 0; // кількість нулів
                    ind = i;
                    while (t <= count_zero_limit && ind < powLen && pow_bin[ind] == '0')
                    {
                        t++;
                        ind++;
                    }

                    if (t > count_zero_limit || ind == powLen)
                        t = -1;

                    if (t >= 0)
                    {
                        for (j = 0; j <= t; j++)
                            res = res * res % mod;
                        res = res * table[t + 1] % mod;
                        i = i + t + 1;
                    }
                    else
                        res = res * table[0] % mod;
                }
                else
                    i++;
            }
            return res;
        }

        public static BigInteger WindowLRMod_Final(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            int i, j, t, powLen;
            Stopwatch stw;
            BigInteger res, temp, found2;
            string pow_bin;
            BigInteger[] table1, table2;

            found = found % mod;
            stw = new Stopwatch();
            stw.Start();
            table1 = new BigInteger[w]; // нульовий елемент саме число, далі 11, 111, ...
            table2 = new BigInteger[w-1]; // нульовий елемент саме число, далі 101, 1001, ...
            temp = table1[0] = table2[0] = found;
            found2 = found;

            BigInteger sqr_found = found * found % mod;
            found = sqr_found * temp % mod;

            table1[1] = found; // записати 11 у двіковій системі

            if (w > 2)
            {
                table2[1] = table1[1] * sqr_found % mod; // 101
                table1[2] = table2[1] * sqr_found % mod; // 111
            }

            for (i = 4; i <= w; i++)
            {
                table2[i-2] = table1[i-2] * sqr_found % mod;

                table1[i-1] = table1[i - 2] * table1[i - 2] % mod;
                table1[i-1] = table1[i-1] * temp % mod;
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            pow_bin = ConvToBinary(pow);
            res = 1;
            powLen = pow_bin.Length;
            i = 0;
            int ind;
            int count_zero_limit = w - 2;
            while (i < powLen)
            {
                if (pow_bin[i] == '1')
                {
                    i++;
                    if ((i == powLen - 1) || (i < powLen - 1 && pow_bin[i] == '0'))
                    {
                        t = 1; // кількість нулів
                        ind = i + t;
                        while (t <= count_zero_limit && ind < powLen && pow_bin[ind] == '0')
                        {
                            t++;
                            ind++;
                        }

                        if (t > count_zero_limit || ind == powLen)
                            t = 0;

                        if (t > 0)
                        {
                            for (j = 1; j <= t + 2; j++)
                                res = res * res % mod;
                            res = res * table2[t] % mod;
                            i = i + t + 1;
                        }
                        else
                        {
                            res = res * res % mod;
                            res = res * table2[0] % mod;
                        }
                    }
                    else
                    {
                        t = 1;  // кількість одиниць
                        ind = i + t;
                        while (t < w && ind < powLen && pow_bin[ind] == '1')
                        {
                            t++;
                            ind++;
                        }

                        for (j = 1; j <= t; j++)
                            res = res * res % mod;
                        res = res * table1[t - 1] % mod;
                        i = i + t - 1;
                    }
                }
                else
                {
                    res = res * res % mod;
                    i++;
                }

            }
            return res;
        }

        public static BigInteger Bonus(BigInteger found, BigInteger po, BigInteger mod, int w, out double table_time)
        {
            found = found % mod;
            BigInteger res = found;//change found
            Stopwatch stw = new Stopwatch();
            stw.Start();


            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            return res;
        }

        public static BigInteger SlideRL(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {

            // 5.1
            found = found % mod;

            int i;
            int count_elem = 1 << (w - 1);
            MyList<BigInteger> table = new MyList<BigInteger>();
            for (i = 1; i <= count_elem; i++)
                table.Add(1);

            BigInteger a = found;
            BigInteger res = 1;

            string pow_bin = ConvToBinary(pow);

            int pow_length = pow_bin.Length;
            i = pow_length - 1;
            int j, index;
            while (i >= 0)
            {
                j = w;
                if (pow_bin[i] == '1')
                {
                    while (i - j + 1 < 0)
                        j--;

                    while (pow_bin[i - j + 1] == '0' && j > 1)
                        j--;

                    i = i - j;
                    string subs = pow_bin.Substring(i + 1, j);
                    index = (Convert.ToInt32(subs, 2) - 1) / 2;
                    if (index >= 0)
                        table[index] = table[index] * a % mod;

                    for (int k = 0; k < j; k++)
                        a = a * a % mod;
                }
                else
                {
                    a = a * a % mod;
                    i--;
                }
            }

            Stopwatch stw = new Stopwatch();
            stw.Start();
            BigInteger temp;
            for (i = count_elem - 1; i > 0; i--)
            {
                table[i - 1] = table[i - 1] * table[i] % mod;
                temp = table[i] * table[i] % mod;
                table[0] = table[0] * temp % mod;
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            return table[0];

            /*

            // залишити закоментованим
            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            BigInteger power = 1 << (w - 1);
            var table = SlideRLTable(found, mod, power, w);

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            BigInteger res = 1;
            BigInteger temp = found;

            string binary = bif.ToBin(pow);

            while (binary.Length > 0)
            {
                int index = binary.Length - 1;
                if (binary.Length < w || binary[index - w + 1] == '0')
                {
                    if (binary[index] == '1')
                        res = mul(res, temp, mod);

                    for (int j = 0; j < table.Count; j++)
                        table[j] = mul(table[j], table[j], mod);

                    temp = mul(temp, temp, mod);

                    binary = binary.Remove(index, 1);
                }
                else
                {
                    int c = Convert.ToInt32(binary.Substring(index - w + 1, w), 2);
                    res = mul(res, table[c - power], mod);

                    temp = mul(temp, table[table.Count - 1], mod);

                    for (int k = 0; k < w; k++)
                        for (int j = 0; j < table.Count; j++)
                            table[j] = mul(table[j], table[j], mod);

                    binary = binary.Remove(index - w + 1, w);
                }
            }
            return res;
            */
        }
        
        public static BigInteger SlideRL_Dic(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {

            // 5.2. Dic
            found = found % mod;

            int i;
            int count_elem = 1 << (w - 1);
            Dictionary<string, BigInteger> table = new Dictionary<string, BigInteger>();
            string temp;
            BigInteger a = found;
            int count_temp = (count_elem << 1) - 1;
            for (i = 1; i <= count_temp; i = i + 2)
            {
                temp = ConvToBinary(i);
                table.Add(temp, 1);
            }


            string pow_bin = ConvToBinary(pow);

            int pow_length = pow_bin.Length;
            i = pow_length - 1;
            int j;

            while (i >= 0)
            {
                j = w;
                if (pow_bin[i] == '1')
                {
                    while (i - j + 1 < 0)
                        j--;

                    while (pow_bin[i - j + 1] == '0' && j > 1)
                        j--;

                    i = i - j;
                    string subs = pow_bin.Substring(i + 1, j);
                    table[subs] = table[subs] * a % mod;

                    for (int k = 0; k < j; k++)
                        a = a * a % mod;
                }
                else
                {
                    a = a * a % mod;
                    i--;
                }
            }



            Stopwatch stw = new Stopwatch();
            stw.Start();
            string temp_previous;
            string str_one = "1";
            BigInteger sqr;
            for (i = count_temp; i > 1; i = i - 2)
            {
                temp = ConvToBinary(i);
                temp_previous = ConvToBinary(i - 2);
                table[temp_previous] = table[temp_previous] * table[temp] % mod;
                sqr = table[temp] * table[temp] % mod;
                table[str_one] = table[str_one] * sqr % mod;
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            return table[str_one];

        }

        public static BigInteger SlideLR(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            MyList<BigInteger> table = new MyList<BigInteger>(); ;

            string temp;
            BigInteger temp_value = found;
            int temp_size = (1 << w) - 1;

            BigInteger sqr_found = found * found % mod;
            int i;
            for (i = 3; i <= temp_size; i = i + 2)
            {
                temp = ConvToBinary(i);
                temp_value = temp_value * sqr_found % mod;
                table.Add(temp_value);
            }

            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            BigInteger res = 1;

            string pow_bin = ConvToBinary(pow);

            i = 0;
            int j, index;
            int pow_length = pow_bin.Length;
            while (i < pow_length)
            {
                j = w;
                if (pow_bin[i] == '1')
                {
                    while (i + j <= pow_length && pow_bin[i + j - 1] == '0' && j > 1)
                        j--;

                    if (j > 1 && i + j <= pow_length)
                    {
                        for (int k = 0; k < j; k++)
                            res = res * res % mod;

                        string subs = pow_bin.Substring(i, j);
                        index = (Convert.ToInt32(subs, 2) - 1) / 2 - 1;
                        res = res * table[index] % mod;
                        i = i + j;
                    }
                    else
                    {
                        res = res * res % mod;
                        res = res * found % mod;
                        i++;
                    }
                }
                else
                {
                    res = res * res % mod;
                    i++;
                }
            }

            return res;

        }

        public static BigInteger SlideLR_Dic(BigInteger found, BigInteger pow, BigInteger mod, int w, out double table_time)
        {
            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            Dictionary<string, BigInteger> table = new Dictionary<string, BigInteger>();

            string temp;
            BigInteger temp_value = found;
            int temp_size = (1 << w) - 1;

            BigInteger sqr_found = found * found % mod;
            int i;
            for (i = 3; i <= temp_size; i = i + 2)
            {
                temp = ConvToBinary(i);
                temp_value = temp_value * sqr_found % mod;
                table.Add(temp, temp_value);
            }

            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            BigInteger res = 1;

            string pow_bin = ConvToBinary(pow);

            i = 0;
            int j;
            int pow_length = pow_bin.Length;
            while (i < pow_length)
            {
                j = w;
                if (pow_bin[i] == '1')
                {
                    while (i + j <= pow_length && pow_bin[i + j - 1] == '0' && j > 1)
                        j = j - 1;
                    if (j > 1 && i + j <= pow_length)
                    {
                        for (int k = 0; k < j; k++)
                            res = res * res % mod;

                        string subs = pow_bin.Substring(i, j);
                        res = res * table[subs] % mod;
                        i = i + j;
                    }
                    else
                    {
                        res = res * res % mod;
                        res = res * found % mod;
                        i++;
                    }
                }
                else
                {
                    res = res * res % mod;
                    i++ ;
                }
            }

            return res;
        }

        /*public static BigInteger SlideLR(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {
            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            BigInteger pow = 2 * ((1 << w) - (BigInteger)Pow((-1), w)) / 3 - 1;
            var table = NAFLRTable(found, mod, pow, w);
            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            MyList<BigInteger> x = power.ToNAF();

            BigInteger res = 1;
            for (BigInteger i = x.Count - 1; i > -1;)
            {
                var max = new MyList<BigInteger>();
                if (x[i] == 0)
                {
                    max.Add(0);
                    max.Add(1);
                }
                else
                    max = FindLargest2(x, i, w);

                for (int d = 0; d < max[1]; d++)
                    res = mul(res, res, mod);

                if (max[0] > 0)
                    res = mul(res, table[(bif.Abs(max[0]) / 2)], mod);
                else if (max[0] < 0)
                    res = mul(res, Euclid_2_1(mod, table[(bif.Abs(max[0]) / 2)]), mod);

                i = i - max[1];
            }
            return res;
        }*/

        public static BigInteger NAFSlideRL(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {
            bool mont = GetMontgomery;

            BigInteger res = 1;
            MyList<int> x = ToNAF(power);

            found = found % mod;

            int sign;
            if (w % 2 == 1)
                sign = -1;
            else
                sign = 1;

            int last_table_element = (((1 << w) - sign) << 1) / 3 - 1;

            MyList<BigInteger> table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();

            int temp;
            BigInteger a = found;

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref a, ref res, mod);

                table.Add(res);
                table_inv.Add(res);

                for (int i = 3; i <= last_table_element; i += 2)
                {
                    table.Add(res);
                    table_inv.Add(res);
                }

                for (int i = 0; i < x.Count;)
                {
                    List<int> max = new List<int>();
                    if (x[i] == 0)
                    {
                        max.Add(0);
                        max.Add(1);
                    }
                    else
                        max = FindLargest1(x, i, w);

                    if (max[0] > 0)
                    {
                        temp = max[0] >> 1;
                        table[temp] = MontgomeryMultDomain(table[temp], a, mod, parameters);
                    }
                    else if (max[0] < 0)
                    {
                        temp = (-max[0]) >> 1;
                        table_inv[temp] = MontgomeryMultDomain(table_inv[temp], a, mod, parameters);
                    }

                    for (int d = 0; d < max[1]; d++)
                        a = MontgomeryMultDomain(a, a, mod, parameters);

                    i = i + max[1];
                }

                Stopwatch stw = new Stopwatch();
                int count_elem = table.Count;
                BigInteger temp_val;

                stw.Start();
                for (int i = count_elem - 1; i > 0; i--)
                {
                    table[i - 1] = MontgomeryMultDomain(table[i - 1], table[i], mod, parameters);
                    temp_val = MontgomeryMultDomain(table[i], table[i], mod, parameters);
                    table[0] = MontgomeryMultDomain(table[0], temp_val, mod, parameters);

                    table_inv[i - 1] = MontgomeryMultDomain(table_inv[i - 1], table_inv[i], mod, parameters);
                    temp_val = MontgomeryMultDomain(table_inv[i], table_inv[i], mod, parameters);
                    table_inv[0] = MontgomeryMultDomain(table_inv[0], temp_val, mod, parameters);
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                res = MontgomeryMultDomain(table[0], MontgomeryInverse(mod, table_inv[0], parameters), mod, parameters);

                res = outMontgomeryDomain(res, mod, parameters);

                return res;
            }
            else
            {
                table.Add(1);
                table_inv.Add(1);

                for (int i = 3; i <= last_table_element; i += 2)
                {
                    table.Add(1);
                    table_inv.Add(1);
                }

                for (int i = 0; i < x.Count;)
                {
                    List<int> max = new List<int>();
                    if (x[i] == 0)
                    {
                        max.Add(0);
                        max.Add(1);
                    }
                    else
                        max = FindLargest1(x, i, w);

                    if (max[0] > 0)
                    {
                        temp = max[0] >> 1;
                        table[temp] = table[temp] * a % mod;
                    }
                    else if (max[0] < 0)
                    {
                        temp = (-max[0]) >> 1;
                        table_inv[temp] = table_inv[temp] * a % mod;
                    }

                    for (int d = 0; d < max[1]; d++)
                        a = a * a % mod;

                    i = i + max[1];
                }


                Stopwatch stw = new Stopwatch();
                int count_elem = table.Count;
                BigInteger temp_val;

                stw.Start();
                for (int i = count_elem - 1; i > 0; i--)
                {
                    table[i - 1] = table[i - 1] * table[i] % mod;
                    temp_val = table[i] * table[i] % mod;
                    table[0] = table[0] * temp_val % mod;

                    table_inv[i - 1] = table_inv[i - 1] * table_inv[i] % mod;
                    temp_val = table_inv[i] * table_inv[i] % mod;
                    table_inv[0] = table_inv[0] * temp_val % mod;
                }
                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                return table[0] * Euclid_2_1(mod, table_inv[0]) % mod;
            }


            /*
            BigInteger res = 1;
            MyList<int> x = power.ToNAF();

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            int sign;
            if (w % 2 == 1)
                sign = -1;
            else
                sign = 1;

            int last_table_element = (((1 << w) - sign) << 1) / 3 - 1;

            MyList<BigInteger> table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();

            table.Add(found);
            table_inv.Add(Euclid_2_1(mod, found));

            BigInteger sqr_found = found * found % mod;

            for (int i = 3; i <= last_table_element; i += 2)
            {
                table.Add(table[(i >> 1) - 1] * sqr_found % mod);
                table_inv.Add(Euclid_2_1(mod, table[i >> 1]));
            }

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;


            for (int i = 0; i < x.Count;)
            {
                List<int> max = new List<int>();
                if (x[i] == 0)
                {
                    max.Add(0);
                    max.Add(1);
                }
                else
                    max = FindLargest1(x, i, w);

                if (max[0] > 0)
                    res = res * table[max[0] >> 1] % mod;
                else if (max[0] < 0)
                    res = res * table_inv[(-max[0]) >> 1] % mod;

                for (int d = 0; d < max[1]; d++)
                    for (int j = 0; j < table.Count; j++)
                    {
                        table[j] = table[j] * table[j] % mod;
                        table_inv[j] = table_inv[j] * table_inv[j] % mod;
                    }

                i = i + max[1];
            }


            return res;
            */


            /*
            BigInteger res = 1;
            MyList<int> x = power.ToNAF();

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            BigInteger pow = 2 * ((1 << w) - (int)Pow((-1), w)) / 3 - 1;
            MyList<BigInteger> table = NAFRLTable(found, mod, pow, w);

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            for (int i = 0; i < x.Count;)
            {
                List<int> max = FindLargest1(x, i, w);

                if (max[0] > 0)
                    res = mul(res, table[(bif.Abs(max[0]) / 2)], mod);
                else if (max[0] < 0)
                    res = mul(res, Euclid_2_1(mod, table[(bif.Abs(max[0]) / 2)]), mod);

                for (int d = 0; d < max[1]; d++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = mul(table[j], table[j], mod);

                i = i + max[1];
            }
            return res;
            */
        }

        public static BigInteger NAFSlideLR(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {
            bool mont = GetMontgomery;

            BigInteger res = 1;
            MyList<int> x = ToNAF(power);

            found = found % mod;

            MyList<BigInteger> table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();
            Stopwatch stw = new Stopwatch();

            if (mont)
            {
                MyList <BigInteger> parameters = MontgomeryMethods.toMontgomeryDomain(ref found, ref res, mod);

                stw.Start();

                int sign;
                if (w % 2 == 1)
                    sign = -1;
                else
                    sign = 1;

                int last_table_element = (((1 << w) - sign) << 1) / 3 - 1;

                table.Add(found);
                table_inv.Add(MontgomeryInverse(mod, found, parameters));

                BigInteger sqr_found = MontgomeryMultDomain(found, found, mod, parameters);

                for (int i = 3; i <= last_table_element; i += 2)
                {
                    table.Add(MontgomeryMultDomain(table[(i >> 1) - 1], sqr_found, mod, parameters));
                    table_inv.Add(MontgomeryInverse(mod, table[i >> 1], parameters));
                }

                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                for (int i = x.Count; i > 0;)
                {
                    List<int> max = new List<int>();
                    if (x[i - 1] == 0 || i == 1)
                    {
                        max.Add(x[i - 1]);
                        max.Add(1);
                    }
                    else
                        max = FindLargest2(x, i - 1, w);

                    for (int d = 0; d < max[1]; d++)
                        res = MontgomeryMultDomain(res, res, mod, parameters);

                    if (max[0] > 0)
                        res = MontgomeryMultDomain(res, table[max[0] >> 1], mod, parameters);
                    else if (max[0] < 0)
                        res = MontgomeryMultDomain(res, table_inv[(-max[0]) >> 1], mod, parameters);

                    i = i - (int) max[1];
                }

                res = outMontgomeryDomain(res, mod, parameters);
            }
            else
            {
                stw.Start();

                int sign;
                if (w % 2 == 1)
                    sign = -1;
                else
                    sign = 1;

                int last_table_element = (((1 << w) - sign) << 1) / 3 - 1;

                table.Add(found);
                table_inv.Add(Euclid_2_1(mod, found));

                BigInteger sqr_found = found * found % mod;

                for (int i = 3; i <= last_table_element; i += 2)
                {
                    table.Add(table[(i >> 1) - 1] * sqr_found % mod);
                    table_inv.Add(Euclid_2_1(mod, table[i >> 1]));
                }

                stw.Stop();
                table_time = stw.Elapsed.TotalMilliseconds;

                for (int i = x.Count; i > 0;)
                {
                    List<int> max = new List<int>();
                    if (x[i - 1] == 0 || i == 1)
                    {
                        max.Add(x[i - 1]);
                        max.Add(1);
                    }
                    else
                        max = FindLargest2(x, i - 1, w);

                    for (int d = 0; d < max[1]; d++)
                        res = res * res % mod;

                    if (max[0] > 0)
                        res = res * table[max[0] >> 1] % mod;
                    else if (max[0] < 0)
                        res = res * table_inv[(-max[0]) >> 1] % mod;

                    i = i - (int) max[1];
                }
            }

            return res;

            /*
                     BigInteger res = 1;
                     MyList<BigInteger> x = power.ToNAF();

                     Stopwatch stw = new Stopwatch();
                     found = found % mod;
                     stw.Start();

                     BigInteger pow = 2 * ((int)Pow(2, w) - (int)Pow((-1), w)) / 3 - 1;
                     MyList<BigInteger> table = NAFLRTable(found, mod, pow, w);

                     stw.Stop();
                     table_time = stw.Elapsed.TotalMilliseconds;

                     for (int i = x.Count - 1; i > -1;)
                     {
                         List<BigInteger> max = new List<BigInteger>();
                         if (x[i] == 0)
                         {
                             max.Add(0);
                             max.Add(1);
                         }
                         else
                             max = FindLargest2(x, i, w);

                         for (int d = 0; d < max[1]; d++)
                             res = mul(res, res, mod);

                         if (max[0] > 0)
                             res = mul(res, table[(bif.Abs(max[0]) / 2)], mod);
                         else if (max[0] < 0)
                             res = mul(res, Euclid_2_1(mod, table[(bif.Abs(max[0]) / 2)]), mod);

                         i = i - (int)max[1];
                     }
                     return res;
             */

        }

        public static BigInteger NAFWindowRL(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {

            BigInteger res = 1;
            MyList<int> x = ToWNAF(power, w);
            MyList<BigInteger> table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();
            int pow = 1 << (w - 1);

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            int count = pow >> 1;
            for (BigInteger i = 0; i < count; i++)
            {
                table.Add(1);
                table_inv.Add(1);
            }

            int temp;
            BigInteger a = found;
            for (int i = x.Count - 1; i > -1; i--)
            {
                if (x[i] > 0)
                {
                    temp = x[i] >> 1;
                    table[temp] = table[temp] * a % mod;
                }
                else if (x[i] < 0)
                {
                    temp = (-x[i]) >> 1;
                    table_inv[temp] = table_inv[temp] * a % mod;
                }

                a = a * a % mod;
            }

            int count_elem = table.Count;
            BigInteger temp_val;
            for (int i = count_elem - 1; i > 0; i--)
            {
                table[i - 1] = table[i - 1] * table[i] % mod;
                temp_val = table[i] * table[i] % mod;
                table[0] = table[0] * temp_val % mod;

                table_inv[i - 1] = table_inv[i - 1] * table_inv[i] % mod;
                temp_val = table_inv[i] * table_inv[i] % mod;
                table_inv[0] = table_inv[0] * temp_val % mod;
            }
            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            return table[0] * Euclid_2_1(mod, table_inv[0]) % mod;


            /*
            BigInteger res = 1;
            MyList<int> x = ToWNAF(power, w);
            MyList<BigInteger> table = new MyList<BigInteger>();
            BigInteger pow = 1 << ( w - 1);

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            for (BigInteger i = 1; i < pow; i += 2)
                table.Add(BinaryRL(found, i, mod));
            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            for (int i = x.Count - 1; i > -1; i--)
            {
                if (x[i] > 0)
                    res = mul(res,  table[(int)(x[i] / 2)], mod);
                else if (x[i] < 0)
                    res = mul(res,  Euclid_2_1(mod, table[(-x[i] / 2)]), mod);

                for (int j = 0; j < table.Count; j++)
                    table[j] = mul(table[j], table[j], mod);
            }
            return res;
            */

        }

        public static BigInteger NAFWindowLR(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {

            BigInteger res = 1;
            var x = ToWNAF(power, w);
            var table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();
            int pow = 1 << (w - 1);

            found = found % mod;
            Stopwatch stw = new Stopwatch();
            stw.Start();
            table.Add(found);
            table_inv.Add(Euclid_2_1(mod, found));

            BigInteger sqr_found = found * found % mod;
            int count = pow >> 1;
            for (BigInteger i = 1; i < count; i++)
            {
                table.Add(table[i - 1] * sqr_found % mod);
                table_inv.Add(Euclid_2_1(mod, table[i]));
            }
            stw.Stop();

            int temp;
            table_time = stw.Elapsed.TotalMilliseconds;
            for (int i = 0; i < x.Count; i++)
            {
                res = res * res % mod;
                
                if (x[i] > 0)
                {
                    temp = x[i] >> 1;
                    res = res * table[temp] % mod;
                }
                else if (x[i] < 0)
                {
                    temp = (-x[i]) >> 1;
                    res = res * table_inv[temp] % mod;
                }
            }
            return res;


            /*
            BigInteger res = 1;
            var x = ToWNAF(power, w);
            var table = new MyList<BigInteger>();
            BigInteger pow = (BigInteger)Pow(2, w - 1);

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();
            for (BigInteger i = 1; i < pow; i += 2)
                table.Add(BinaryLR(found, i, mod));
            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            for (BigInteger i = 0; i < x.Count; i++)
            {
                res = mul(res, res, mod);
                if (x[i] > 0)
                    res = mul(res, table[(x[i] / 2)], mod);
                else if (x[i] < 0)
                    res = mul(res, Euclid_2_1(mod, table[(-x[i] / 2)]), mod);
            }
            return res;
            */

        }

        public static BigInteger wNAFSlideRL(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {
            BigInteger res = 1;
            var x = ToWNAF(power, w);
            found = found % mod;

            var table = new MyList<BigInteger>();
            BigInteger pow = 2 * ((BigInteger)Pow(2, w) - (BigInteger)Pow((-1), w)) / 3 - 1;
            Stopwatch stw = new Stopwatch();
            stw.Start();

            for (BigInteger i = 1; i <= pow; i += 2)
                table.Add(BinaryRL(found, i, mod));

            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;

            for (int i = x.Count - 1; i > -1; )
            {
                var max = FindLargest2(x, i, w);

                if (max[0] > 0)
                    res = mul(res, table[(BigInteger)(Abs((decimal)x[i]) / 2)], mod);
                else if (max[0] < 0)
                    res = mul(res, Euclid_2_1(mod, table[(BigInteger)(Abs((decimal)x[i]) / 2)]), mod);

                for (int d = 0; d < max[1]; d++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = mul(table[j], table[j], mod);

                i = i - max[1];
            }
            return res;
        }

        public static BigInteger wNAFSlideLR(BigInteger found, BigInteger power, BigInteger mod, int w, out double table_time)
        {

            BigInteger res = 1;
            var x = ToWNAF(power, w);

            Stopwatch stw = new Stopwatch();
            found = found % mod;
            stw.Start();

            int sign;
            if (w % 2 == 1)
                sign = -1;
            else
                sign = 1;

            int last_table_element = (((1 << w) - sign) << 1) / 3 - 1;

            MyList<BigInteger> table = new MyList<BigInteger>();
            MyList<BigInteger> table_inv = new MyList<BigInteger>();

            table.Add(found);
            table_inv.Add(Euclid_2_1(mod, found));

            BigInteger sqr_found = found * found % mod;

            for (int i = 3; i <= last_table_element; i += 2)
            {
                table.Add(table[(i >> 1) - 1] * sqr_found % mod);
                table_inv.Add(Euclid_2_1(mod, table[i >> 1]));
            }

            stw.Stop();
            table_time = stw.Elapsed.TotalMilliseconds;

            for (int i = 0; i < x.Count;)
            {
                List<int> max = new List<int>();
                if (x[i] == 0 || i == 0)
                {
                    max.Add(x[i]);
                    max.Add(1);
                }
                else
                    max = FindLargest1(x, i, w);

                for (int d = 0; d < max[1]; d++)
                    res = res * res % mod;

                if (max[0] > 0)
                    res = res * table[max[0] >> 1] % mod;
                else if (max[0] < 0)
                    res = res * table_inv[(-max[0]) >> 1] % mod;

                i = i + max[1];
            }
            return res;


            /*
            BigInteger res = 1;
            var x = ToWNAF(power,w);
            found = found % mod;

            var table = new MyList<BigInteger>();
            BigInteger pow = 2 * ((BigInteger)Pow(2, w) - (BigInteger)Pow((-1), w)) / 3 - 1;
            Stopwatch stw = new Stopwatch();
            stw.Start();
            for (BigInteger i = 1; i <= pow; i += 2)
                table.Add(BinaryLR(found, i, mod));
            stw.Stop();

            table_time = stw.Elapsed.TotalMilliseconds;
            for (int i = 0; i < x.Count; )
            {
                var max = new MyList<int>();
                if (x[i] == 0)
                {
                    max.Add(0);
                    max.Add(1);
                }
                else
                    max = FindLargest2(x, i, w);

                for (int d = 0; d < max[1]; d++)
                    res = mul(res, res, mod);

                if (max[0] > 0)
                    res = mul(res, table[(BigInteger)(Abs((decimal)max[0]) / 2)], mod);
                else if (max[0] < 0)
                    res = mul(res, Euclid_2_1(mod, table[(BigInteger)(Abs((decimal)max[0]) / 2)]), mod);

                i = i + max[1];
            }
            return res;
            */

        }

        #endregion
    }
    abstract class GenFunctions
    {
        public static List<int> ReadString(string s)
        {
            string st = "";
            List<string> sts = new List<string>();
            for (int i = 0; i < s.Length; i++)
            {
                if (s[i] == ' ')
                    continue;
                if (s[i] != ';')
                    st = st + s[i];
                else
                {
                    sts.Add(st);
                    st = "";
                }
            }
            if (st != "")
            {
                sts.Add(st);
                st = "";
            }

            List<int> x = new List<int>();
            List<int> pows = new List<int>();
            for (int i = 0; i < sts.Count; i++)
            {
                x.Clear();
                for (int j = 0; j < sts[i].Length; j++)
                {
                    if (sts[i][j] != '-')
                        st = st + sts[i][j];
                    if (sts[i][j] == '-')
                    {
                        x.Add(Convert.ToInt16(st));
                        st = "";
                    }
                }
                if (st != "")
                {
                    x.Add(Convert.ToInt16(st));
                    st = "";
                }

                if (x.Count == 1)
                    pows.Add(x[0]);
                else
                {
                    int step;
                    if (x.Count == 2)
                        step = 1;
                    else
                        step = x[2];
                    for (int j = x[0]; j <= x[1]; j += step)
                        pows.Add(j);
                }
            }
            pows.Sort();
            for (int i = 0; i < pows.Count - 1; i++)
            {
                if (pows[i] == pows[i + 1])
                {
                    pows.RemoveAt(i);
                    i--;
                }
            }
            return pows;
        }
        public static int rand(int l, int h)
        {
            Random r = new Random();
            int x = r.Next(l, h);
            return x;
        }
        public delegate BigInteger random_num(int n, string type);

        public static BigInteger random_max(int n, string type = "Bytes")
        {
            var rng = new RNGCryptoServiceProvider();
            byte[] bytes = new byte[n + 1];
            rng.GetBytes(bytes);
            if (type == "Bytes")
            {
                bytes[n - 1] |= 0x80;//128
                bytes[n] = 0;
                return new BigInteger(bytes);
            }
            else
            {
                string s = "01";
                for (int i = 0; i < n - 1; i++)
                {
                    s += bytes[i] % 2 == 0 ? '1' : '0';
                }

                return new BigInteger(Convert.ToInt32(s, 2));
            }
        }

        public static BigInteger random_max_bites(int n)//n - length in bit
        {
            string s = "";
            var rand = new Random();
            for (int i = 0; i < n + 1; i++)
            {
                s += rand.Next(2) == 1 ? '1' : '0';
            }

            BigInteger b = new BigInteger(Convert.ToInt32(s, 2));

            return b;
        }

        public static BigInteger random_two(int n, string type)//n - length in bytes
        {
            if (type == "Bytes")
            {
                var rng = new RNGCryptoServiceProvider();
                byte[] bytes = new byte[n + 1];
                bytes[n] = 1;
                for (int i = 0; i < bytes.Length - 1; i++)
                {
                    bytes[i] = 0;
                }

                return new BigInteger(bytes);
            }
            else
            {
                string s = "";
                s += '1';
                for (int i = 1; i < n + 1; i++)
                {
                    s += '0';
                }
                return new BigInteger(Convert.ToInt32(s, 2));
            }

        }
        public static BigInteger random_simple(int n, string type)//n - length in bytes
        {
            BigInteger b;
            again:
            b = random_odd(n, type);
            if (type == "Bytes")
            {
                if (b.IsProbablePrime())
                    return b;
            }
            else
            {
                if (b.IsProbablePrime())
                    //if (Primality_Tests.Ferma(b))
                    return b;
            }
            goto again;
        }
        public static BigInteger random_odd(int n, string type)//n - length in bytes
        {
            var rng = new RNGCryptoServiceProvider();
            byte[] bytes = new byte[n + 1];
            rng.GetBytes(bytes);
            if (type == "Bytes")
            {
                bytes[n - 1] |= 0x80;
                bytes[n] = 0;
                bytes[0] |= 1;
                if (bytes[0] % 2 == 0)
                {
                    bytes[0]++;
                }

                return new BigInteger(bytes);
            }
            else
            {
                char[] mas = new char[n + 1];
                mas[n - 1] = '1';
                mas[n] = '0';
                mas[0] = '1';
                for (int i = 1; i < n - 1; i++)
                {
                    mas[i] += bytes[i] % 2 == 0 ? '1' : '0';
                }
                Array.Reverse(mas);
                return new BigInteger(Convert.ToInt32(new string(mas), 2));
            }

        }
    }

    public static class BigIntegerExtensions
    {
        //Тест Миллера-Рабина на простоту
        public static bool IsProbablePrime(this BigInteger source, int certainty = 5)
        {
            if (source == 2 || source == 3)
                return true;
            if (source < 2 || source % 2 == 0)
                return false;

            BigInteger d = source - 1;
            int s = 0;

            while (d % 2 == 0)
            {
                d /= 2;
                s += 1;
            }

            // There is no built-in method for generating random BigInteger values.
            // Instead, random BigIntegers are constructed from randomly generated
            // byte arrays of the same length as the source.
            RandomNumberGenerator rng = RandomNumberGenerator.Create();
            byte[] bytes = new byte[source.ToByteArray().LongLength];
            BigInteger a;

            for (int i = 0; i < certainty; i++)
            {
                do
                {
                    rng.GetBytes(bytes);
                    a = new BigInteger(bytes);
                }
                while (a < 2 || a >= source - 2);

                BigInteger x = BigInteger.ModPow(a, d, source);
                if (x == 1 || x == source - 1)
                    continue;

                for (int r = 1; r < s; r++)
                {
                    x = BigInteger.ModPow(x, 2, source);
                    if (x == 1)
                        return false;
                    if (x == source - 1)
                        break;
                }

                if (x != source - 1)
                    return false;
            }

            return true;
        }
        public static BigInteger Pow(this BigInteger value, BigInteger exponent)
        {
            BigInteger _base = 1;
            for (BigInteger i = 0; i < exponent; i++)
            {
                _base *= value;
            }
            return _base;
        }
        public static BigInteger Euclid_2_1(this BigInteger mod, BigInteger found)
        {
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
            return inv;
        }
        public static MyList<int> ToNAF(this BigInteger k)
        {
            MyList<int> mas_k = new MyList<int>();
            int i = 0;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    mas_k.Add((int)(2 - (k % 4)));
                    k = k - mas_k[i];
                }
                else
                    mas_k.Add(0);


                k = k / 2;
                i++;
            }
            return mas_k;
        }
        public static BigInteger TwoPow(this BigInteger pow)
        {
            BigInteger result = 1;
            for (BigInteger i = 0; i < pow; i++)
            {
                result = result << 1;
            }
            return result;
        }
        private static string ConvToBinary(this BigInteger num)
        {
            string x = "";
            while (num != 0)
            {
                x += num % 2;
                num /= 2;
            }
            char[] charArray = x.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }
        //Тест Миллера-Рабина на простоту код студента
        public static bool Prime_Test_Miller_Rabin(this BigInteger y)
        {
            BigInteger b = 0, T, a;
            bool flag = true;
            int k = 0, i = 0; //a - osnova
            if ((y & 1) == 0) return false;  // nuzhno ne4etnoe 4islo
            if (y == 3) return true;
            do
            {
                int N = (int)Ceiling(Ceiling(BigInteger.Log(y - 1, 2)) / 8); //kol-vo byte v 4isle (y-1)
                int rnd = GenFunctions.rand(1, N);                //
                a = GenFunctions.random_max(rnd);          //generiruem osnovu <<a>> v diapazone ot 1 do (y-1)
            }
            while ((a >= y - 1) || (a <= 1));


            if (BigInteger.GreatestCommonDivisor(a, y) != 1) return false;
            else
            {
                while (flag == true)
                {
                    k++;
                    b = (y - 1) >> k;                    // y-1 = 2^k * b,
                    if ((b & 1) == 1) flag = false;      // b - ne4etnoe
                }
                T = BigInteger.ModPow(a, b, y);
                if ((T == 1) || (T == y - 1)) return true;

                else
                {
                    for (i = 1; i < k; i++)
                    {
                        T = BigInteger.ModPow(T, 2, y);
                        if (T == y - 1) return true;
                    }
                }
            }

            return false;
        }
    }

}

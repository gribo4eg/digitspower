using System;
using System.Collections.Generic;
using System.Numerics;
using System.Security.Cryptography;
using System.Diagnostics;

namespace Polinoms
{
    class MyList<T> : List<T>
    {
        //public MyList() : base() { }
        public T this[BigInteger index]    // Indexer declaration
        {
            get { return base[(int)index]; }
            set { base[(int)index] = value; }
        }
    }
    abstract class Pow
    {
       
        private static MyList<BigInteger> Table(BigInteger found, BigInteger pow, int w, BigInteger mod) 
        {
            var table = new MyList<BigInteger>();

            table.Add(found);
            for (BigInteger i = 0; i < BigInteger.Parse((Math.Pow(2, w) - 2).ToString()); i++)
                table.Add((table[i] * found) % mod);//Приведення до степеня після кожного кроку
            return table;
        }

        private static MyList<BigInteger> SlideRLTable(BigInteger found, BigInteger pow, int w, BigInteger mod) 
        {
            var table = new MyList<BigInteger>();
            table.Add(BinaryRL(found,pow,mod));
            for (int i = 0; i < pow - 1; i++)
                table.Add(table[i] * found % mod);
            return table;
        }
        private static MyList<BigInteger> SlideLRTable(BigInteger found, BigInteger pow, int w, BigInteger mod) 
        {
            var table = new MyList<BigInteger>();
            table.Add(BinaryLR(found,pow,mod));
            for (int i = 0; i < pow - 1; i++)
                table.Add(table[i] * found % mod);
            return table;
        }
        private static MyList<BigInteger> NAFRLTable(BigInteger found, BigInteger pow, int w, BigInteger mod) 
        {
            var table = new MyList<BigInteger>();

            for (int i = 1; i <= pow; i += 2)
                table.Add(BinaryRL(found,i,mod));
            return table;
        }
        private static MyList<BigInteger> NAFLRTable(BigInteger found, BigInteger pow, int w, BigInteger mod) 
        {
            var table = new MyList<BigInteger>();

            for (int i = 1; i <= pow; i += 2)
                table.Add(BinaryLR(found,i,mod));
            return table;
        }
         public static BigInteger BinaryRL(BigInteger found, BigInteger pow, BigInteger mod)
        {
            BigInteger res = 1;
            BigInteger t = found;

            string Binary = ConvToBinary(pow);
            for (int i = Binary.Length - 1; i >= 0; i--)
            {
                if (Binary[i] == '1')
                    res = (t * res) % mod;//Приведення до степеня після кожного кроку
                t = (t * t) % mod;
            }
            return res;
        }

        public static cpolinom BinaryRL(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom t = cpolinom.Copy(p);
            string Binary = Convert.ToString(n, 2);
            for (int i = Binary.Length - 1; i >= 0; i--)
            {
                if (Binary[i] == '1')
                    res = t * res % ir;
                t = t * t % ir;
            }
            return res;
        }
        public static cpolinom BinaryLR(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom t = cpolinom.Copy(p);
            string Binary = Convert.ToString(n, 2);
            for (int i = 0; i < Binary.Length; i++)
            {
                res = res * res % ir;
                if (Binary[i] == '1')
                    res = t * res % ir;
            }
            return res;
        }

        private static List<cpolinom> Table(cpolinom p, cpolinom ir, int w)
        {
            List<cpolinom> table = new List<cpolinom>();
            table.Add(p % ir);
            for (int i = 0; i < Math.Pow(2, w) - 2; i++)
                table.Add(table[i] * p % ir);
            return table;
        }


        private static List<string> windows(string s, int w)
        {
            string st = s;

            while (st.Length % w != 0)
                st = "0" + st;


            List<string> bins = new List<string>();
            for (int i = 0; i < st.Length / w; i++)
                bins.Add(st.Substring(i * w, w));

            return bins;
        }

        public static cpolinom WindowRL(cpolinom p, cpolinom ir, int n, int w)
        {
            List<cpolinom> table = Table(p, ir, w);
            cpolinom res = new cpolinom("1", p.mod);

            List<string> bins = windows(Convert.ToString(n, 2), w);

            for (int i = bins.Count - 1; i > -1; i--)
            {
                int c = Convert.ToInt32(bins[i], 2);
                if (c != 0) res = res * table[c - 1] % ir;

                for (int k = 0; k < w; k++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = table[j] * table[j] % ir;
            }
            return res;
        }
        public static cpolinom WindowLR(cpolinom p, cpolinom ir, int n, int w)
        {
            List<cpolinom> table = Table(p, ir, w);
            cpolinom res = new cpolinom("1", p.mod);

            List<string> bins = windows(Convert.ToString(n, 2), w);

            for (int i = 0; i < bins.Count; i++)
            {
                for (int k = 0; k < w; k++)
                    res = res * res % ir;

                int c = Convert.ToInt32(bins[i], 2);
                if (c != 0) res = res * table[c - 1] % ir;
            }
            return res;
        }

        public static BigInteger SlideRL(BigInteger found, BigInteger n, BigInteger mod, int w, out float table_time)
        {
            int power = (int)Math.Pow(2, w - 1);

            Stopwatch stw = new Stopwatch();
            stw.Start();

            var table = SlideRLTable(found, power, w, mod);

            stw.Stop();
            table_time = (float)(stw.ElapsedTicks / 10);

            BigInteger res = 1;
            BigInteger temp = found;
            //List<string> bins = windows(ConvToBinary(n), w);
            string binary = ConvToBinary(n);

            while (binary.Length > 0)
            {
                int index = binary.Length - 1;
                if (binary.Length < w || binary[index - w + 1] == '0')
                {
                    if (binary[index] == '1')
                        res = (res * temp) % mod;

                    for (int j = 0; j < table.Count; j++)
                        table[j] = table[j] * table[j] % mod;

                    temp = temp * temp % mod;

                    binary = binary.Remove(index, 1);
                }
                else
                {
                    int c = Convert.ToInt32(binary.Substring(index - w + 1, w), 2);
                    res = res * table[c - power] % mod;

                    temp = temp * table[table.Count - 1] % mod;

                    for (int k = 0; k < w; k++)
                        for (int j = 0; j < table.Count; j++)
                            table[j] = table[j] * table[j] % mod;

                    binary = binary.Remove(index - w + 1, w);
                }
            }
            return res;
        }
        public static BigInteger SlideLR(BigInteger found, BigInteger n, BigInteger mod, int w, out float table_time)
        {
            int power = (int)Math.Pow(2, w - 1);

            Stopwatch stw = new Stopwatch();
            stw.Start();
            var table = SlideLRTable(found, power, w ,mod);
            stw.Stop();
            table_time = (float)(stw.ElapsedTicks / 10);

            BigInteger res = 1;
            string binary = ConvToBinary(n);

            while (binary.Length > 0)
            {
                if (binary.Length < w || binary[0] == '0')
                {
                    res = res * res % mod;
                    if (binary[0] == '1')
                        res = res * found % mod;
                    binary = binary.Remove(0, 1);
                }
                else
                {
                    int c = Convert.ToInt32(binary.Substring(0, w), 2);

                    for (int k = 0; k < w; k++)
                        res = res * res % mod;

                    res = res * table[c - power] % mod;
                    binary = binary.Remove(0, w);
                }
            }
            return res;
        }

        private static List<int> ToNAF(int n)
        {
            List<int> res = new List<int>();
            int k = n;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    res.Add(2 - k % 4);
                    k = k - res[res.Count - 1];
                }
                else
                    res.Add(0);
                k = k / 2;
            }
            return res;
        }

        public static cpolinom NAFBinaryRL(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom c = cpolinom.Copy(p);
            List<int> x = ToNAF(n);
            for (int i = 0; i < x.Count; i++)
            {
                if (x[i] == 1) res = res * c % ir;
                else if (x[i] == -1) res = res * c.EvkPolinom(ir) % ir;
                c = c * c % ir;
            }
            return res;
        }
        public static cpolinom NAFBinaryLR(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom c = cpolinom.Copy(p);
            List<int> x = ToNAF(n);
            for (int i = x.Count - 1; i > -1; i--)
            {
                res = res * res % ir;
                if (x[i] == 1) res = res * c % ir;
                else if (x[i] == -1) res = res * c.EvkPolinom(ir) % ir;
            }
            return res;
        }
        public static cpolinom meth7_2(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom c = cpolinom.Copy(p);
            int k = n;
            while (k >= 1)
            {
                int temp = 0;
                if (k % 2 != 0)
                {
                    temp = 2 - k % 4;
                    k -= temp;
                }

                if (temp == 1) res = res * c % ir;
                else if (temp == -1) res = res * c.EvkPolinom(ir) % ir;

                k = k / 2;
                c = c * c % ir;

            }
            return res;
        }


        private static List<int> FindLargest1(List<int> x, int i, int w)
        {
            int j = i;
            int pow = 1;
            int temp = x[i];

            int max_j = i;
            int max = x[i];
            while (j - i + 1 < w && j < x.Count - 1)
            {
                j++;

                pow = pow * 2;
                temp = temp + pow * x[j];

                if (temp % 2 != 0)
                {
                    max_j = j;
                    max = temp;
                }

            }
            max_j = max_j - i + 1;

            List<int> r = new List<int>();
            r.Add(max);
            r.Add(max_j);
            return r;
        }
        private static List<int> FindLargest2(List<int> x, int i, int w)
        {
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

            List<int> r = new List<int>();
            r.Add(max);
            r.Add(max_j);
            return r;
        }

        public static cpolinom NAFSlideRL(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToNAF(n);
            int pow = 2 * ((int)Math.Pow(2, w) - (int)Math.Pow((-1), w)) / 3 - 1;
            List<cpolinom> table = NAFRLTable(p, ir, pow, w);

            for (int i = 0; i < x.Count; )
            {
                List<int> max = FindLargest1(x, i, w);

                if (max[0] > 0)
                    res = res * table[(int)(Math.Abs(max[0]) / 2)] % ir;
                else if (max[0] < 0)
                    res = res * table[(int)(Math.Abs(max[0]) / 2)].EvkPolinom(ir) % ir;

                for (int d = 0; d < max[1]; d++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = table[j] * table[j] % ir;

                i = i + max[1];
            }
            return res;
        }
        public static cpolinom NAFSlideLR(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToNAF(n);
            int pow = 2 * ((int)Math.Pow(2, w) - (int)Math.Pow((-1), w)) / 3 - 1;
            List<cpolinom> table = NAFLRTable(p, ir, pow, w);

            for (int i = x.Count - 1; i > -1; )
            {
                List<int> max = new List<int>();
                if (x[i] == 0)
                {
                    max.Add(0);
                    max.Add(1);
                }
                else
                    max = FindLargest2(x, i, w);

                for (int d = 0; d < max[1]; d++)
                    res = res * res % ir;

                if (max[0] > 0)
                    res = res * table[(int)(Math.Abs(max[0]) / 2)] % ir;
                else if (max[0] < 0)
                    res = res * table[(int)(Math.Abs(max[0]) / 2)].EvkPolinom(ir) % ir;

                i = i - max[1];
            }
            return res;
        }

        public static cpolinom NAFWindowRL(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToNAF(n);
            int pow = (int)Math.Pow(2, w - 1);
            List<cpolinom> table = NAFRLTable(p, ir, pow, w);

            for (int i = 1; i < pow; i += 2)
                table.Add(BinaryRL(p, ir, i));

            for (int i = x.Count - 1; i > -1; i--)
            {
                if (x[i] > 0)
                    res = res * table[(int)(x[i] / 2)] % ir;
                else if (x[i] < 0)
                    res = res * table[(int)(-x[i] / 2)].EvkPolinom(ir) % ir;

                for (int j = 0; j < table.Count; j++)
                    table[j] = table[j] * table[j] % ir;
            }
            return res;
        }
        public static cpolinom NAFWindowLR(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToNAF(n);
            int pow = (int)Math.Pow(2, w - 1);
            List<cpolinom> table = NAFLRTable(p, ir, pow, w);

            for (int i = 0; i < x.Count; i++)
            {
                res = res * res % ir;
                if (x[i] > 0)
                    res = res * table[(int)(x[i] / 2)] % ir;
                else if (x[i] < 0)
                    res = res * table[(int)(-x[i] / 2)].EvkPolinom(ir) % ir;
            }
            return res;
        }

        private static List<int> ToWNAF(int n, int w)
        {
            List<int> res = new List<int>();
            int k = n;
            int r;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    r = k % (int)Math.Pow(2, w);

                    if (r >= Math.Pow(2, w - 1))
                        res.Insert(0, r - (int)Math.Pow(2, w));
                    else
                        res.Insert(0, r);

                    k = k - res[0];
                }
                else
                    res.Insert(0, 0);
                k = k / 2;
            }
            return res;
        }

        public static cpolinom wNAFSlideRL(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToWNAF(n, w);
            int pow = 2 * ((int)Math.Pow(2, w) - (int)Math.Pow((-1), w)) / 3 - 1;

            List<cpolinom> table = NAFRLTable(p, ir, pow, w);

            for (int i = x.Count - 1; i > -1; )
            {
                List<int> max = FindLargest1(x, i, w);

                if (max[0] > 0)
                    res = res * table[(int)(Math.Abs(x[i]) / 2)] % ir;
                else if (max[0] < 0)
                    res = res * table[(int)(Math.Abs(x[i]) / 2)].EvkPolinom(ir) % ir;

                for (int d = 0; d < max[1]; d++)
                    for (int j = 0; j < table.Count; j++)
                        table[j] = table[j] * table[j] % ir;

                i = i - max[1];
            }
            return res;
        }
        public static cpolinom wNAFSlideLR(cpolinom p, cpolinom ir, int n, int w)
        {
            cpolinom res = new cpolinom("1", p.mod);
            List<int> x = ToWNAF(n, w);
            int pow = 2 * ((int)Math.Pow(2, w) - (int)Math.Pow((-1), w)) / 3 - 1;

            List<cpolinom> table = NAFLRTable(p, ir, pow, w);

            for (int i = 0; i < x.Count; )
            {
                List<int> max = new List<int>();
                if (x[i] == 0)
                {
                    max.Add(0);
                    max.Add(1);
                }
                else
                    max = FindLargest2(x, i, w);

                for (int d = 0; d < max[1]; d++)
                    res = res * res % ir;

                if (max[0] > 0)
                    res = res * table[(int)(Math.Abs(x[i]) / 2)] % ir;
                else if (max[0] < 0)
                    res = res * table[(int)(Math.Abs(x[i]) / 2)].EvkPolinom(ir) % ir;

                i = i + max[1];
            }
            return res;
        }

        public static cpolinom AddSubRL(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom c = cpolinom.Copy(p);

            string BinaryT = Convert.ToString(3 * n, 2);
            string BinaryN = Convert.ToString(n, 2);
            while (BinaryN.Length < BinaryT.Length) BinaryN = "0" + BinaryN;

            for (int i = BinaryN.Length - 2; i >= 0; i--)
            {
                if (BinaryT[i] != '0' && BinaryN[i] == '0') res = res * c % ir;
                else if (BinaryT[i] == '0' && BinaryN[i] != '0') res = res * c.EvkPolinom(ir) % ir;
                c = c * c % ir;
            }

            return res;
        }
        public static cpolinom AddSubLR(cpolinom p, cpolinom ir, int n)
        {
            cpolinom res = new cpolinom("1", p.mod);
            cpolinom c = cpolinom.Copy(p);

            string BinaryT = Convert.ToString(3 * n, 2);
            string BinaryN = Convert.ToString(n, 2);
            while (BinaryN.Length < BinaryT.Length) BinaryN = "0" + BinaryN;

            for (int i = 0; i < BinaryN.Length - 1; i++)
            {
                res = res * res % ir;
                if (BinaryT[i] != '0' && BinaryN[i] == '0') res = res * c % ir;
                else if (BinaryT[i] == '0' && BinaryN[i] != '0') res = res * c.EvkPolinom(ir) % ir;
            }

            return res;
        }



        private static double LogDifBase(int Base, int argument)
        {
            return Math.Log(argument) / Math.Log(Base);
        }
        private static int get_number_bit(int value, int Base)
        {
            int number_bit = (int)Math.Floor(LogDifBase(Base, value));
            if (LogDifBase(Base, value) > number_bit || number_bit == 0)
                number_bit++;
            return number_bit;
        }
        private static int[] ReCompute(int a, int b, int a_max, int b_max, int k)
        {
            int na = a;
            int nb = b;
            int temp;
            if (a > a_max)
            {
                temp = get_number_bit((int)Math.Pow(2, (na - a_max)), 3) - 1;
                nb = nb + temp;
                if (nb > b_max) nb = b_max;
                if (a_max > 0) temp = a_max - 1;
                else temp = 0;

                if ((Math.Abs(k - Math.Pow(2, a_max) * Math.Pow(3, nb)) < Math.Abs(k - Math.Pow(2, temp) * Math.Pow(3, nb + 1))) || nb == b_max)
                    na = a_max;
                else
                {
                    na = temp;
                    b++;
                }
            }
            int[] x = new int[2];
            x[0] = na;
            x[1] = nb;
            return x;

        }
        private static int[] BestApproximation1(int k, int a_max, int b_max)
        {
            int min_x = a_max;
            double y = Math.Round(-min_x * LogDifBase(3, 2) + LogDifBase(3, k));
            if (y > b_max)
                y = b_max;
            else if (y < 0)
                y = 0;
            double mindelta = Math.Abs(y + min_x * LogDifBase(3, 2) - LogDifBase(3, k));
            int x = 0;
            for (x = 0; x <= a_max; x++)
            {
                y = Math.Round(-x * LogDifBase(3, 2) + LogDifBase(3, k));
                if (y > b_max)
                    y = b_max;
                else if (y < 0)
                    y = 0;

                double delta = Math.Abs(y + x * LogDifBase(3, 2) - LogDifBase(3, k));
                if (mindelta > delta)
                {
                    min_x = x;
                    mindelta = delta;
                }
            }
            int a = min_x;
            int b = (int)Math.Round(-min_x * LogDifBase(3, 2) + LogDifBase(3, k));
            if (b > b_max)
                b = b_max;
            int[] r = new int[2];
            r[0] = a;
            r[1] = b;
            return r;
        }
        private static List<int[]> ToDBNS1RL(int n, int a_max, int b_max)
        {
            List<int[]> t = new List<int[]>();
            int i = 0;
            int s = 1;
            int k = n;
            int am;
            int bm;
            while (k > 0)
            {
                t.Add(new int[3]);
                int[] x = BestApproximation1(k, a_max, b_max);
                am = x[0];
                bm = x[1];
                int z = (int)Math.Pow(2, am) + (int)Math.Pow(3, bm);
                t[i][0] = s;
                t[i][1] = am;
                t[i][2] = bm;

                if (k < z)
                    s = -s;
                k = Math.Abs(k - z);

                i++;
            }
            return t;
        }
        public static cpolinom DBNS1RL(cpolinom p, cpolinom ir, int n)
        {
            int a_max = 15;
            int b_max = 17;
            List<int[]> x = ToDBNS1RL(n, a_max, b_max);
            cpolinom t = cpolinom.Copy(p);
            cpolinom res = new cpolinom("1", p.mod);

            for (int j = 0; j < x[0][2]; j++)
                t = t * t * t;
            for (int j = 0; j < x[0][1]; j++)
                t = t * t;
            if (x[0][0] == -1)
                res = t.EvkPolinom(ir);
            else if (x[0][0] == 1)
                res = t;
            for (int i = 0; i < x.Count - 1; i++)
            {
                int u = x[i][1] - x[i + 1][1];
                int v = x[i][2] - x[i + 1][2];
                for (int j = 0; j <= u; j++)
                    t = t * t;
                for (int j = 0; j <= v; j++)
                    t = t * t * t;
                if (x[0][0] == -1)
                    res = res * t.EvkPolinom(ir);
                else if (x[0][0] == 1)
                    res = res * t;
            }
            return res;
        }
        //19
        //20

        private static List<int[]> ToDBNS2RL(int n)
        {
            List<int[]> x = new List<int[]>();
            int v = n;
            int i = 0;
            while (v != 0)
            {
                x.Insert(0, new int[3]);
                int a = 0;
                while (v % 2 == 0 && v != 0)
                {
                    a++;
                    v /= 2;
                }
                int b = 0;
                while (v % 3 == 0 && v != 0)
                {
                    b++;
                    v /= 3;
                }
                x[0][1] = a;
                x[0][2] = b;

                v--;
                if (v != 0)
                {
                    if (v % 6 != 0)
                    {
                        v += 2;
                        x[0][0] = -1;
                    }
                    else
                        x[0][0] = 1;
                }
                else
                    x[0][0] = 1;
                i++;
            }
            return x;
        }
        public static cpolinom DBNS2RL(cpolinom p, cpolinom ir, int n)
        {
            List<int[]> x = ToDBNS2RL(n);

            cpolinom res = new cpolinom("1", p.mod);
            cpolinom t = cpolinom.Copy(p);
            for (int i = x.Count - 1; i > -1; i--)
            {
                for (int j = 0; j < x[i][1]; j++)
                    t = t * t % ir;
                for (int j = 0; j < x[i][2]; j++)
                    t = t * t * t % ir;

                if (x[i][0] == 1)
                    res = res * t % ir;
                else if (x[i][0] == -1)
                    res = res * t.EvkPolinom(ir) % ir;
            }
            return res;
        }
        private static List<int[]> ToDBNS2LR(int n)
        {
            List<int[]> x = new List<int[]>();
            x.Add(new int[3]);
            int i = 0;
            int v = n;
            x[i][0] = 0;
            while (v != 0)
            {
                int a = 0;
                while (v % 2 == 0 && v != 0)
                {
                    a++;
                    v = v / 2;
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
        public static cpolinom DBNS2LR(cpolinom p, cpolinom ir, int n)
        {
            List<int[]> x = ToDBNS2LR(n);
            cpolinom res = cpolinom.Copy(p);

            for (int i = 1; i < x.Count; i++)
            {
                for (int j = 0; j < x[i][1]; j++)
                    res = res * res % ir;
                for (int j = 0; j < x[i][2]; j++)
                    res = res * res * res % ir;

                if (x[i][0] == 1)
                    res = res * p % ir;
                else if (x[i][0] == -1)
                    res = res * p.EvkPolinom(ir) % ir;
            }
            for (int j = 0; j < x[0][1]; j++)
                res = res * res % ir;
            for (int j = 0; j < x[0][2]; j++)
                res = res * res * res % ir;
            return res;
        }
    }

}

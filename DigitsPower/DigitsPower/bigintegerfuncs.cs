using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace DigitsPower
{
    class bif
    {
        public static BigInteger Abs(BigInteger dec)
        {
            if (dec < 0)
                return -dec;
            else
                return dec;
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
        public static List<BigInteger> wNAF(BigInteger n, int w)
        {
            List<BigInteger> res = new List<BigInteger>();
            BigInteger k = n;
            BigInteger r;
            while (k >= 1)
            {
                if (k % 2 != 0)
                {
                    r = k % (int)Math.Pow(2, w);

                    if (r >= (int)Math.Pow(2, w - 1))
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
            return (BigInteger.Log(_base)/Math.Log(argument));
        }
    }
}

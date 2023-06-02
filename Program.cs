using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using static System.Net.Mime.MediaTypeNames;

namespace PullingPoZavodu
{
    internal class Program
    {
        static void Main(string[] args)
        {
            double outerDiameter = 10;
            double outerDiameterTolerance = 0.5;
            double innerDiameter = 5;
            double innerDiameterTolerance = 0.25;
            double height = 2;
            double tapeThickness = 0.1;
            double topAdditionalThickness = 0.05;
            double bottomAdditionalThickness = 0.05;

            double d0KolpNesht = Functions.d0KolpNesht(outerDiameter, outerDiameterTolerance, innerDiameter, innerDiameterTolerance, height, tapeThickness, topAdditionalThickness, bottomAdditionalThickness);

            double d0KolpSht = Functions.d0KolpSht(outerDiameter, outerDiameterTolerance, innerDiameter, innerDiameterTolerance, height, tapeThickness, topAdditionalThickness, bottomAdditionalThickness, 0.5);

            Console.WriteLine("The diameter of the non-threaded cylindrical object is {0}", d0KolpNesht);
            Console.WriteLine("The diameter of the threaded cylindrical object is {0}", d0KolpSht);
            Console.ReadKey();
        }
    }



    public class Functions
    {
        public static double d0KolpNesht(double outerDiameter, double outerDiameterTolerance, double innerDiameter, double innerDiameterTolerance, double height, double tapeThickness, double topAdditionalThickness, double bottomAdditionalThickness)
        {
            double dn = outerDiameter + outerDiameterTolerance / 2;
            double dw = innerDiameter + innerDiameterTolerance / 2;
            double h = height;
            double s = tapeThickness + (topAdditionalThickness + bottomAdditionalThickness) / 2;
            double pi = Math.PI;
            double v = pi * (dn * dn - dw * dw) * (h - s) / 4 + pi * dw * dw * s / 4;
            double d0KolpNesht = Math.Sqrt(4 * v / pi / s);
            d0KolpNesht = Math.Round(d0KolpNesht, 2);
            return d0KolpNesht;
        }

        public static double d0KolpSht(double outerDiameter, double outerDiameterTolerance, double innerDiameter, double innerDiameterTolerance, double height, double tapeThickness, double topAdditionalThickness, double bottomAdditionalThickness, double radius)
        {
            double dn = outerDiameter + outerDiameterTolerance / 2;
            double dw = innerDiameter + innerDiameterTolerance / 2;
            double h = height;
            double s = tapeThickness + (topAdditionalThickness + bottomAdditionalThickness) / 2;
            double r = radius;
            double pi = Math.PI;
            double vcil = pi * (dn * dn - dw * dw) * (h - s - r) / 4;
            double vdno = pi * s * (dw - 2 * r) * (dw - 2 * r) / 4;
            double vtor = pi * pi * (dw / 2 - r) * ((r + s) * (r + s) - r * r) / 2;
            double v = vcil + vdno + vtor;
            double d0KolpSht = Math.Sqrt(4 * v / pi / s);
            d0KolpSht = Math.Round(d0KolpSht, 2);
            return d0KolpSht;
        }

        public static void KoeffWyt()
        {
            // Get the input values
            double md0 = float.Parse(Console.ReadLine());
            double d0 = float.Parse(Console.ReadLine());
            double ms0 = float.Parse(Console.ReadLine());
            double sl = float.Parse(Console.ReadLine());
            double delsw = float.Parse(Console.ReadLine());
            double delsn = float.Parse(Console.ReadLine());
            double dn = float.Parse(Console.ReadLine());
            double dw = float.Parse(Console.ReadLine());
            double deldn = float.Parse(Console.ReadLine());
            double deldw = float.Parse(Console.ReadLine());
            double s0 = sl + delsw / 2 + delsn / 2;

            // Calculate the initial values
            double m1 = float.Parse(Console.ReadLine());
            int n = 0;
            if (md0 - m1 >= 0)
            {
                m1 = md0;
                double mu1 = ms0;
                n = 1;
            }
            else
            {
                // Get the input coefficient of thinning
                double mu1 = float.Parse(Console.ReadLine());

                // Get the input average coefficient of extraction
                double mdsr = float.Parse(Console.ReadLine());

                // Calculate the number of extractions
                double nsr = (Math.Log(md0) - Math.Log(m1)) / Math.Log(mdsr) + 1;
                n = (int)Math.Round(nsr);
                if ((nsr - n) > 0.3)
                {
                    n = n + 1;
                }
            }
            Console.WriteLine("Количество вытяжек: " + (n - 1));

            // Create arrays to store the draw ratios, average diameters, thinning ratios, wall thicknesses, matrix diameters, and punch diameters.
            double[] md = new double[n];
            double[] dsr = new double[n];
            double[] ms = new double[n];
            double[] s = new double[n];
            double[] dm = new double[n];
            double[] dp = new double[n];
            // Initialize the arrays.
            for (int i = 0; i < n; i++)
            {
                md[i] = m1;
                ms[i] = ms0;
                s[i] = s0;
            }

            // Calculate the draw ratios, average diameters, thinning ratios, wall thicknesses, matrix diameters, and punch diameters for each draw.
            for (int i = 1; i <= n; i++)
            {
                if (i == 1)
                {
                    md[i - 1] = m1;
                    ms[i - 1] = ms0;
                }
                else
                {
                    md[i - 1] = md[i - 2] * m1;
                    ms[i - 1] = ms[i - 2] * ms0;
                }
                dsr[i - 1] = (md[i - 1] + 1) / 2;
                s[i - 1] = s[i - 2] * ms[i - 1];
                dm[i - 1] = dn + (dsr[i - 1] - 1) * deldn;
                dp[i - 1] = dw - (dsr[i - 1] - 1) * deldw;
            }

            // Display the results.
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("Вытяжка " + (i + 1) + ":");
                Console.WriteLine("  md = " + md[i]);
                Console.WriteLine("  dsr = " + dsr[i]);
                Console.WriteLine("  ms = " + ms[i]);
                Console.WriteLine("  s = " + s[i]);
                Console.WriteLine("  dm = " + dm[i]);
                Console.WriteLine("  dp = " + dp[i]);
            }
        }

        //расчет диам. матриц(кроме 1-й и окончат.) по отчету ЛВМИ 1960
        public static double DmSmirAl(double dmf, double dpf, double dps, double stdef)
        {
            double a = Math.Sqrt(3) / 2;
            double b = 1 - 1 / Math.Sqrt(3);
            double c = a * (stdef - b * Math.Log(dpf / dps));
            double d = Math.Exp(c);
            return Math.Round(Math.Sqrt(dps * dps + (dmf * dmf - dpf * dpf) / d), 2);
        }

        //суммар. степень деформ по Смирнову-Аляеву
        public static double StepDefSum(double d0, double s0, double dn, double dw)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(d0 / dw) + b * Math.Log(4 * d0 * s0 / (dn * dn - dw * dw));
        }

        //степ. деформ. на свертке по Смирнову-Аляеву
        public static double StepDefSwert(double d0, double s0, double dn, double dw)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(d0 / dw) + b * Math.Log(4 * d0 * s0 / (dn * dn - dw * dw));
        }

        //степ. деформ. на вытяжке по Смирнову-Аляеву
        public static double StepDefWyt(double dmf, double dms, double dpf, double dps)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(dpf / dps) + b * Math.Log((dmf * dmf - dpf * dpf) / (dms * dms - dps * dps));
        }

        public static double Hh(double dz, double s, double dp, double dm, double r)
        {
            double Pi = Math.PI;
            return (dz * dz * s - (dp - 2 * r) * s * (Pi * (2 * r + s) + dp - 2 * r)) / (dm * dm - dp * dp) + s + r;
        }

        public static double Hn(double dz, double s, double dp, double dm)
        {
            return s * ((dz * dz - dp * dp) / (dm * dm - dp * dp) + 1);
        }
    }

    public class Sub
    {
        public static void DiamKru(string[] args)
        {
            // Get the input values
            int outerDiameter = int.Parse(Console.ReadLine());
            int innerDiameter = int.Parse(Console.ReadLine());
            int height = int.Parse(Console.ReadLine());
            int thickness = int.Parse(Console.ReadLine());
            int topMargin = int.Parse(Console.ReadLine());
            int bottomMargin = int.Parse(Console.ReadLine());

            // Calculate the required diameter of the circle
            int circleDiameter = CalculateCircleDiameter(outerDiameter, innerDiameter, height, thickness, topMargin, bottomMargin);

            // Print the output
            Console.WriteLine("The required diameter of the circle is: " + circleDiameter);
        }

        private static int CalculateCircleDiameter(int outerDiameter, int innerDiameter, int height, int thickness, int topMargin, int bottomMargin)
        {
            // Calculate the minimum and maximum heights of the detail
            int minHeight = CalculateMinimumHeight(outerDiameter, innerDiameter, thickness, topMargin, bottomMargin);
            int maxHeight = CalculateMaximumHeight(outerDiameter, innerDiameter, thickness, topMargin, bottomMargin);

            // Calculate the required diameter of the circle
            int circleDiameter = 0;
            if (minHeight == maxHeight)
            {
                circleDiameter = CalculateCircleDiameterForSingleHeight(outerDiameter, innerDiameter, minHeight);
            }
            else
            {
                circleDiameter = CalculateCircleDiameterForMultipleHeights(outerDiameter, innerDiameter, minHeight, maxHeight);
            }

            return circleDiameter;
        }

        private static int CalculateMinimumHeight(int outerDiameter, int innerDiameter, int thickness, int topMargin, int bottomMargin)
        {
            // Calculate the minimum height of the detail
            int minHeight = (outerDiameter - innerDiameter) / 2 + topMargin + bottomMargin;

            return minHeight;
        }

        private static int CalculateMaximumHeight(int outerDiameter, int innerDiameter, int thickness, int topMargin, int bottomMargin)
        {
            // Calculate the maximum height of the detail
            int maxHeight = (outerDiameter - innerDiameter) / 2 + thickness + topMargin + bottomMargin;

            return maxHeight;
        }

        private static int CalculateCircleDiameterForSingleHeight(int outerDiameter, int innerDiameter, int height)
        {
            // Calculate the required diameter of the circle for a single height
            int circleDiameter = (outerDiameter + innerDiameter) / 2 * height / outerDiameter;

            return circleDiameter;
        }

        private static int CalculateCircleDiameterForMultipleHeights(int outerDiameter, int innerDiameter, int minHeight, int maxHeight)
        {
            // Calculate the required diameter of the circle for multiple heights
            int circleDiameter = (outerDiameter + innerDiameter) / 2 * (minHeight + maxHeight) / outerDiameter;

            return circleDiameter;
        }
    }


    public class SmirAl
    {
        public static void RadPuans()
        {
            // Get the radius and number of steps from the user.
            double r = double.Parse(Console.ReadLine());
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the intermediate diameters and the radii.
            double[] md = new double[n];
            double[] dsr = new double[n];
            double[] rp = new double[n];

            // Populate the arrays with the values from the user.
            for (int i = 0; i < n; i++)
            {
                dsr[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the radii for each step.
            for (int k = 0; k < n; k++)
            {
                if (k == n - 1)
                {
                    rp[k] = r;
                }
                else
                {
                    rp[k] = dsr[k] * (1 - md[k + 1]) / 2;
                    rp[k] = Math.Round(rp[k], 1);
                }

                // Print the value of the radius.
                Console.WriteLine(rp[k]);
            }
        }

        public static void WysotaPoluf()
        {
            // Get the radius and number of steps from the user.
            double r = double.Parse(Console.ReadLine());
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the intermediate diameters, the radii, the coefficients of drawing, the modulus of the matrix and the height.
            double[] dsr = new double[n];
            double[] md = new double[n];
            double[] rp = new double[n];
            double[] ms = new double[n];
            double[] h = new double[n];
            double[] x = new double[n];
            double[] y = new double[n];

            // Populate the arrays with the values from the user.
            for (int i = 0; i < n; i++)
            {
                dsr[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
                ms[i] = double.Parse(Console.ReadLine());
                rp[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the heights for each step.
            for (int k = 0; k < n; k++)
            {
                x[0] = 1;
                y[0] = 1;
                x[k] = x[k - 1] * ms[k];
                y[k] = y[k - 1] * md[k];
                h[k] = 0.25 * dsr[k] * (1 / y[k] * y[k] - 1 - 2.28 * rp[k] / dsr[k] + 0.56 * (rp[k] / dsr[k]) * (rp[k] / dsr[k])) / x[k] + rp[k];
                h[k] = Math.Round(h[k], 1);
                Console.WriteLine(h[k]);
            }
        }

        public static void VysotaOtkusa()
        {
            // Get the input values from the user.
            double dmp = double.Parse(Console.ReadLine());
            double dpp = double.Parse(Console.ReadLine());
            double dmo = double.Parse(Console.ReadLine());
            double dpo = double.Parse(Console.ReadLine());
            double al = double.Parse(Console.ReadLine());
            double d0 = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double Pi = Math.PI;
            double h = double.Parse(Console.ReadLine());

            // Calculate the initial volume of the blank.
            double v0 = Pi * d0 * d0 * s / 4;

            // Calculate the volume of the final part.
            double vd = Pi * ((dmo * dmo - dpo * dpo) * h + dpo * dpo * s) / 4;

            // Calculate the cone angle in radians.
            double alf = Pi * al / 180;

            // Calculate the height of the conical part of the cut.
            double hk = (dmp - dmo) / (2 * Math.Tan(alf));
            hk = Math.Round(hk, 1);

            // Calculate the length of the cylindrical part of the cut.
            double lk = (dmp - dmo) / (2 * Math.Sin(alf));
            lk = Math.Round(lk, 1);

            // Calculate the radius of the punch shoulder.
            double sp = (dmp - dpp) / 2;

            // Calculate the volume of the conical part of the cut.
            double vkon = Pi * lk * (dmp - sp) * (dmp - dpp) / 4 + Pi * lk * (dmo - sp) * (dmp - dpp) / 4;

            // Calculate the remaining volume to be cut.
            double votk = v0 - vd;

            // Check if there is enough metal to cut.
            if (votk < 0)
            {
                Console.WriteLine("нет металла на откус");
                return;
            }

            // Check if the conical part of the cut is large enough.
            if (vkon >= votk)
            {
                // Calculate the length of the cut.
                double ds = dmo - sp;
                double a = Pi * Math.Sin(alf) * sp;
                double b = Pi * ds * sp;
                double c = -votk;
                if (b * b - 4 * a * c < 0)
                {
                    Console.WriteLine("наверно неправильно введены данные");
                    return;
                }
                double d = Math.Sqrt(b * b - 4 * a * c);
                double l = (d - b) / (2 * a);
                l = Math.Round(l, 1);
                Console.WriteLine("длина откуса=" + l);
            }
            else
            {
                // Calculate the height of the cylindrical part of the cut.
                double vc = votk - vkon;
                double hc = 4 * vc / (Pi * (dmp * dmp - dpp * dpp));
                hc = Math.Round(hc, 1);
                Console.WriteLine("высота цил.части откус=" + hc);
                Console.WriteLine("высот. конич.части отк=" + hk);
                double hs = hk + hc;
                Console.WriteLine("суммар высота откуса=" + hs);
            }
        }

        public static void KoefИзвестнИнструм()
        {
            // Get the input values from the user.
            int n = int.Parse(Console.ReadLine());
            double d0 = double.Parse(Console.ReadLine());
            double dn = double.Parse(Console.ReadLine());
            double deln = double.Parse(Console.ReadLine());
            double dw = double.Parse(Console.ReadLine());
            double delw = double.Parse(Console.ReadLine());
            double r = double.Parse(Console.ReadLine());
            double t = double.Parse(Console.ReadLine());
            double tw = double.Parse(Console.ReadLine());
            double tn = double.Parse(Console.ReadLine());

            // Create arrays to store the intermediate diameters, the radii, the coefficients of drawing, the modulus of the matrix and the height.
            double[] dsr = new double[n];
            double[] md = new double[n];
            double[] ms = new double[n];
            double[] s = new double[n];
            double[] dm = new double[n];
            double[] dp = new double[n];
            double[] rp = new double[n];
            double[] h = new double[n];
            double[] x = new double[n];
            double[] y = new double[n];

            // Populate the arrays with the values from the user.
            for (int i = 0; i < n; i++)
            {
                rp[i] = double.Parse(Console.ReadLine());
                dm[i] = double.Parse(Console.ReadLine());
                dp[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the thicknesses of the walls.
            for (int i = 0; i < n; i++)
            {
                s[i] = (dm[i] - dp[i]) / 2;
                if (i == 0 && (dm[i] - dp[i]) / 2 > t)
                {
                    s[i] = (t + tw / 2 + tn / 2);
                }
                dsr[i] = (dm[i] + dp[i]) / 2;
                if (i == 0)
                {
                    md[i] = dsr[i] / d0;
                    ms[i] = s[i] / (t + tw / 2 + tn / 2);
                }
                else
                {
                    md[i] = dsr[i] / dsr[i - 1];
                    ms[i] = s[i] / s[i - 1];
                }
                if (ms[i] > 1)
                {
                    ms[i] = 1;
                }
                Console.WriteLine(md[i]);
                Console.WriteLine(ms[i]);
                Console.WriteLine(dsr[i]);
                Console.WriteLine(s[i]);
            }
            double md0 = dsr[n - 1] / d0;
            double ms0 = s[n - 1] / (t + tw / 2 + tn / 2);
            Console.WriteLine(md0);
            Console.WriteLine(ms0);
            for (int k = 1; k <= n; k++)
            {
                x[0] = 1;
                y[0] = 1;
                x[k] = x[k - 1] * ms[k];
                y[k] = y[k - 1] * md[k];
                h[k] = 0.25 * dsr[k] * (1 / y[k] * y[k] - 1 - 2.28 * rp[k] / dsr[k] + 0.56 * (rp[k] / dsr[k]) * (rp[k] / dsr[k])) / x[k] + rp[k];
            }
            Console.WriteLine("h:");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(h[i]);
            }
            double sd = t / d0;
            Console.WriteLine(sd);
        }

        public static void dpSmirAl()
        {
            // Get the input values
            Console.WriteLine("Введи диам. пуанс. свертки");
            double dpsw = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи диам. пуанс. окончат.выт ");
            double dpow = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи количество вытяжек ");
            int n = int.Parse(Console.ReadLine());

            // Create an array to store the diameters
            double[] dp = new double[n];

            // Set the first and last elements of the array
            dp[0] = dpsw;
            dp[n - 1] = dpow;

            // Calculate the intermediate elements of the array
            for (int i = 1; i < n - 1; i++)
            {
                dp[i] = dp[i - 1] / Math.Pow((dpsw / dpow), (1.0 / n));
                dp[i] = Math.Round(dp[i], 2);
            }

            // Print the array to the console
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(dp[i]);
            }
        }

        public static void одноконМатрПрижим(string[] args)
        {
            // Get the input values
            double pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            double alf = double.Parse(Console.ReadLine());
            double al = pi * alf / 180;
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double mum = 0.05;

            // Calculate the diameter of the cone entrance
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(pi * alf / 180) + 1);
            dk = Math.Round(dk, 1);

            // Calculate the radius of the cone entrance
            double rw = 3 * s;
            double rws = rw + s / 2;

            // Calculate the height of the cone
            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);

            // Calculate the reduction ratio of the diameter
            double md12 = d1 / dk;
            double md11 = dk / d0;

            // Calculate the reduction ratio of the length
            double psi = 1 - dk / d0;

            // Calculate the angle of the frustum
            double fik = pi * (90 - alf) / 180;

            // Calculate the reduction ratio of the cross-section
            double psisr = 1 - Math.Sqrt(md1);

            // Calculate the yield stress in the shell
            double sigms = sigmb * (psisr / psir) * Math.Pow(psir / (1 - psir), 1.0 / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);

            // Calculate the coefficient a
            double a = Math.Log(1 / md11) - psi + s / (2 * rws * Math.Sqrt(md11));

            // Calculate the coefficient b
            double b = 1 + mum * fik;

            // Calculate the coefficient c
            double c = 1 - 18 * sd / (1 - md1);

            // Calculate the yield stress in the ring
            double sigmr = 1.1 * sigms * (a * b / (1 - 0.2 * mum * b * c / md1) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);

            if (ms1 == 1)
            {
                // Calculate the required force for rolling without thinning
                double pbu = Math.PI * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);

                // Set the output value for the required force
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки без утонения=" + pbu + " кг");
            }
            else
            {
                // Calculate the yield stress in the shell after thinning
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);

                // Calculate the yield stress in the ring after thinning
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);

                // Calculate the required force for rolling with thinning
                double pu = Math.PI * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);

                // Set the output value for the required force
                Console.WriteLine("усилие свертки с утонением=" + pu + " кг");

                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;

                // Calculate the yield stress in the shell at the final stage of rolling
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);

                // Calculate the yield stress in the ring at the final stage of rolling
                double sigmzk = sigmazk(mum, msk, al, sigms2k); ;
                sigmzk = Math.Round(sigmzk, 1);

                // Calculate the required force for rolling at the final stage of thinning
                double pk = Math.PI * d1 * s * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmzk = " + sigmzk + " kg / mm^2");

                // Set the output value for the required force
                Console.WriteLine("усилие свертки с утонением в конечн. стадии=" + pk + " кг");
            }

            // Check for buckling
            double ustKonus = (Math.Sqrt((20 * sd) * (20 * sd) * (1 - Math.Sin(al)) + Math.Sin(al) * Math.Sin(al)) - 20 * sd) / Math.Sin(al);
            if (md1 < ustKonus)
            {
                Console.WriteLine("Однокон. матр.с плоск.и КОНИЧ. прижимом");
            }
            else
            {
                Console.WriteLine("Одноконус. матрица с прижимом");
            }

            // Print the input parameters
            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf + " град");
        }

        public static double sigmaz(double mum, double sigmr, double sigms2, double al, double ms1)
        {
            return ((1 + mum * (1 - sigmr / (1.15 * sigms2)) / Math.Sin(al) - mum * Math.Log(1 / ms1) / Math.Sin(al)) * Math.Log(1 / ms1) + sigmr / (1.15 * sigms2) + Math.Sin(al) / 2) * 1.15 * sigms2;
        }

        public static double sigmas2(double sigmb, double md1, double ms1, double psir)
        {
            return sigmb * Math.Pow((1 - md1 * ms1) / psir, psir / (1 - psir));
        }

        public static double sigmazk(double mum, double msk, double al, double sigms2k)
        {
            return 1.15 * sigms2k * ((1 + mum * (1 - Math.Log(1 / msk)) / Math.Sin(al)) * Math.Log(1 / msk) + Math.Sin(al) / 2);
        }

        static void одноконМатрБезПриж()
        {
            // Define the variables
            double pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            double alf = double.Parse(Console.ReadLine());
            double al = pi * alf / 180;
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double mum = 0.05;

            // Calculate the values
            double dk = Math.Round(0.9 * d0, 1);
            double hk = Math.Round((dk - dm1) / 2 / Math.Tan(al));
            double rw = Math.Round(0.05 * d0);
            double rws = rw + s / 2;
            double dkk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(pi * alf / 180) + 1);
            double md12 = d1 / dkk;
            double md11 = dkk / d0;
            double psi = 1 - dkk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * (psisr / psir) * (psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(al)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);

            // Calculate the forces
            double pbu = 0;
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                pbu = pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = pi * d1 * s1 * sigmz;
                pu = Math.Round(pbu, 0);
                //конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
            }

            // Print the results
            Console.WriteLine("Одноконус. матрица без прижима");
            Console.WriteLine("Диаметр входн.кромки конуса={0}", dk);
            Console.WriteLine("Радиус входн.кромки конуса={0}", rw);
            Console.WriteLine("Высота конуса={0}", hk);
            Console.WriteLine("Угол конуса={0}", alf);
            Console.WriteLine("sigms={0} кг / мм^2", sigms);
            Console.WriteLine("sigmr={0} кг / мм^2", sigmr);
            Console.WriteLine("усилие свертки={0} кг", pbu);
        }

        public static void двухконМатр()
        {
            // Define constants
            double mum = 0.05;
            double Pi = Math.PI;

            // Read input values
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate recommended angle for upper cone
            if (sd > 0.012 && sd < 0.018)
            {
                string werhugol = "30 degrees";
                Console.WriteLine("Recommended angle for upper cone = " + werhugol);
            }
            else if (sd > 0.018 && sd < 0.05)
            {
                string werhugol = "45 degrees";
                Console.WriteLine("Recommended angle for upper cone = " + werhugol);
            }

            double alfw = double.Parse(Console.ReadLine());
            double alfn = double.Parse(Console.ReadLine());

            // Convert angles to radians
            double alw = Pi * alfw / 180;
            double alf = Pi * alfn / 180;

            double dk = d1 * Math.Sqrt((1 / (md1 * md1) - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double dw = 0.9 * d0;
            dw = Math.Round(dw, 1);
            double hw = (dw - dk) / (2 * Math.Tan(alw));
            hw = Math.Round(hw, 1);
            double rw = 0.05 * d0;
            rw = Math.Round(rw, 1);
            double rs = (d0 - dm1) / 3;
            rs = Math.Round(rs, 1);
            double r = (d0 - d1) / 2;
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double fik = (alw - alf) / 2;
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (((1 + mum / Math.Tan(alw)) * (Math.Log(1 / md11) - psi) + s / (2 * d1 * Math.Sqrt(md11))) * (1 + mum * fik) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);

            if (ms1 == 1)
            {
                // усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("Force of rolling = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("Final stage");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("Force of rolling = " + pk + " kg");
            }
            Console.WriteLine("sigms = " + sigms + " kg / mm^2");
            Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");

            Console.WriteLine("Двухконус. матрица");
            Console.WriteLine("Диам. верхнего конуса", dw);
            Console.WriteLine("Диаметр нижнего конуса", dk);
            Console.WriteLine("Высота верхнего конуса", hw);
            Console.WriteLine("Высота нижнего конуса", hk);
            Console.WriteLine("Угол верхнего конуса", alfw);
            Console.WriteLine("угол нижнего конуса", alfn);
            Console.WriteLine("Радиус входн.верх конуса", rw);
            Console.WriteLine("Радиус сопряж. конусов", rs);
        }

        public static void радМатрПриж()
        {
            // Define the constants.
            const double mum = 0.05;
            const double Pi = Math.PI;

            // Get the values from the Excel spreadsheet.
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate the radius of the matrix.
            double rms = 0.69 * d1 * (0.16 * Math.Sqrt(18.3 / md1 * md1 + 21.2 + 10.2 * (rps / d1) * (rps / d1) - 41.3 * rps / d1) - 1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);

            // Calculate the reduced modulus of deformation.
            double md11 = (d1 + 2 * rms) / d0;
            double md12 = md1 / md11;
            double alr = Math.Atan2(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)), rm + s1);
            double fi = Pi / 2 - alr;
            double al = alr / 2;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * (psisr / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);

            // Calculate the shear modulus.
            double a = Math.Log(1 / md1) + s / (4 * rm);
            double b = 1 + 1.5 * mum;
            double c = 0.2 * mum * b / md1;
            double d = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * a * b / (1 - c * d);
            sigmr = Math.Round(sigmr, 1);

            // Calculate the rolling force without thinning.
            double pbu;

            // Calculate the rolling force with thinning.
            if (ms1 == 1)
            {
                //усилие при свертке без утонения.
                pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                //конечная стадия
                // Calculate the rolling force with thinning in the final stage.
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }

            Console.WriteLine("Радиальная матрица с прижимом");
            if (sd < (1 - psir * Math.Log(1 / md1) - md1) / 20)
            {
                Console.WriteLine("Требуется плоский и торцевой прижим");
            }
            else
            {
                Console.WriteLine("Плоский прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }

        public static void радМатрБезПриж()
        {
            // Define the variables.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate the radius of the matrix.
            double rms = 0.69 * d1 * (0.16 * Math.Sqrt(18.3 / md1 * md1 + 21.2 + 10.2 * (rps / d1) * (rps / d1) - 41.3 * rps / d1) - 1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);

            // Calculate the other parameters.
            double md11 = (d1 + 2 * rms) / d0;
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double fi = Pi / 2 - alr;
            double al = alr / 2;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * (psisr / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double a = Math.Log(1 / md1) + s / (4 * rm);
            double b = 1 + 1.5 * mum;
            double c = 0.2 * mum * b / md1;
            double d = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * a * b / (1 - c * d);
            sigmr = Math.Round(sigmr, 1);

            // Calculate the required forces.
            double pbu, pu;
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }

            // Print the results.
            Console.WriteLine("Радиальная матрица с прижимом");
            if (sd < (1 - psir * Math.Log(1 / md1) - md1) / 20)
            {
                Console.WriteLine("Требуется плоский и торцевой прижим");
            }
            else
            {
                Console.WriteLine("Плоский прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }

        public static void ВырСвертка(string[] args)
        {

            // Get the values from the user
            Console.WriteLine("Введите тип матрицы: матрица конич.-K, радиальная-R");
            string tipMatriz = Console.ReadLine();
            double sd = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Check if the squeeze is needed
            if (sd > (1 - md1) / 18)
            {
                Console.WriteLine("прижим не нужен на кон. и рад. матр");
            }
            else if (sd > (1 - md1) / 36)
            {
                Console.WriteLine("на конич. матр. прижим не нужен");
            }
            else
            {
                Console.WriteLine("нужен прижим");
            }

            // Check the type of matrix
            if (tipMatriz == "k")
            {
                if (sd > 0.012 && sd < 0.05)
                {
                    Console.WriteLine("Рекомендуется двухконусная матрица");
                }
                else if (sd > 0.05)
                {
                    Console.WriteLine("Рекомендуется одноконусная матрица с большим радиусом");
                }

                // Check if the squeeze is needed for conical matrix
                if (sd <= (1 - md1) / 36)
                {
                    Console.WriteLine("одноконМатрПрижим");
                }
                else
                {
                    Console.WriteLine("матрица однокон.-O, двухкон.-D, однокон. с радиус-OR");
                    string tipKonMatriz = Console.ReadLine();

                    if (tipKonMatriz == "o")
                    {
                        Console.WriteLine("одноконМатрБезПриж");
                    }
                    else if (tipKonMatriz == "d")
                    {
                        Console.WriteLine("двухконМатр");
                    }
                    else if (tipKonMatriz == "or")
                    {
                        Console.WriteLine("одноконМатрРадиус");
                    }
                }
            }
            else if (tipMatriz == "r")
            {
                if (sd > (1 - md1) / 18)
                {
                    Console.WriteLine("радМатрБезПриж");
                }
                else
                {
                    Console.WriteLine("радМатрПрижим");
                }
            }
        }

        public static void одноконМатрРадиус(string[] args)
        {
            // Define constants
            double mum = 0.05;
            double Pi = Math.PI;

            // Get input values
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            //Введи угол  конуса
            double alfn = double.Parse(Console.ReadLine());
            double alf = Pi * alfn / 180;

            // Calculate values
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double md11 = dk / d0;
            double a = md11 + sd;
            double b = (1 - Math.Sin(alf)) * Math.Tan(alf);
            double hm = (a * (1 - b) + b - sd * ms1 - md1) * d0 / (2 * Math.Tan(alf));
            hm = Math.Round(hm, 1);
            double c = hm - hk;
            if (c < 0)
            {
                Console.WriteLine("высота матр. меньше высоты конуса");
                return;
            }
            double dkm = dm1 + 2 * hm * Math.Tan(alf);
            dkm = Math.Round(dkm, 1);
            double md12 = d1 / dk;
            double rw = (d0 - dk - s) / 2;
            rw = Math.Round(rw, 1);
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir))) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(alf)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            if (ms1 == 1)
            {
                // усилие при свертке без утонения
                double pbu = Math.PI * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("Диам.  конуса=" + dkm);
            Console.WriteLine("Высота конуса=" + hm);
            Console.WriteLine("угол  конуса=" + alfn);
            Console.WriteLine("радиус конуса=" + rw);
        }

        public static void УголКонусМатр(string[] args)
        {
            // Define constants
            double Pi = Math.PI;

            // Get input values
            int n = int.Parse(Console.ReadLine());
            double mum = 0.05;

            // Create arrays to store the input values
            double[] s = new double[n - 1];
            double[] d = new double[n - 1];
            double[] md = new double[n - 1];

            // Get the input values and store them in the arrays
            for (int i = 0; i < n - 1; i++)
            {
                s[i] = double.Parse(Console.ReadLine());
                d[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the angles
            for (int i = 0; i < n - 1; i++)
            {
                double a = Math.Sqrt(0.76 * mum * Math.Log(1 / md[i]) * Math.Sqrt(d[i] / s[i]));
                double b = Math.Asin(a);
                double alf = 180 * b / Pi;
                alf = Math.Round(alf, 0);
                Console.WriteLine("alf={0}", alf);
            }
        }
        public static void ВысотаКонусаМатр(string[] args)
        {
            // Define constants
            double Pi = Math.PI;

            // Get input values
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the input values
            double[] dm = new double[n - 1];
            double[] alf = new double[n - 1];

            // Get the input values and store them in the arrays
            for (int i = 0; i < n - 1; i++)
            {
                dm[i] = double.Parse(Console.ReadLine());
            }

            // Get the input radius
            //Введи радиус вход. конуса матр
            double r = double.Parse(Console.ReadLine());

            // Calculate the heights
            for (int i = 1; i < n - 1; i++)
            {
                alf[i] = double.Parse(Console.ReadLine());
                double al = Pi * alf[i] / 180;
                double hm = (dm[i - 1] - dm[i]) / 2 / Math.Tan(al) + r * (1 - Math.Sin(al));
                double dkm = dm[i - 1] + 2 * r * (1 - Math.Sin(al)) * Math.Tan(al);
                dkm = Math.Round(dkm, 1);
                hm = Math.Round(hm, 1);
                Console.WriteLine("hm={0}", hm);
                Console.WriteLine("dkm={0}", dkm);
                Console.WriteLine("r={0}", r);
            }
        }
        public static void ВытБезУтон()
        {
            // Define constants
            double Pi = Math.PI;
            // Calculate the stress
            double psir = 1 - Math.Sqrt(md);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir))) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            Console.WriteLine("sigms={0}", sigms);
        }

        static void УсилиеВытяжки()
        {
            // Define constants
            const double Pi = Math.PI;
            const double mum = 0.05;

            // Get the number of iterations
            int n = int.Parse(Console.ReadLine());

            // Get the material properties
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Create arrays to store the data
            double[] alf = new double[n - 1];
            double[] md = new double[n - 1];
            double[] ms = new double[n - 1];
            double[] s = new double[n];
            double[] dsr = new double[n];

            // Read the data from the console
            for (int i = 0; i < n - 1; i++)
            {
                alf[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
                ms[i] = double.Parse(Console.ReadLine());
            }
            for (int i = 0; i < n; i++)
            {
                s[i] = double.Parse(Console.ReadLine());
                dsr[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the results
            for (int i = 0; i < n - 1; i++)
            {
                // Convert the angle from degrees to radians
                double alfaRad = alf[i] * Pi / 180;

                if (ms[i] < 1)
                {
                    // Calculate the yield stresses
                    double sigms1 = sigmas1(sigmb, md[i], psir);
                    double sigms2 = sigmas2(sigmb, md[i], ms[i], psir);
                    sigms2 = Math.Round(sigms2, 1);

                    // Calculate the reduced yield stresses
                    double sigmr1 = 1.1 * sigms1 * sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);
                    double sigr1 = sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);

                    // Calculate the bursting stresses
                    double sigmz = sigmaz2(ms[i], alfaRad, sigr1, sigms2);
                    sigmz = Math.Round(sigmz, 1);

                    // Calculate the bursting forces
                    double pst = pste(dsr[i], s[i], sigmz);
                    double ptr = Pi * ptre(ms[i], dsr[i], s[i], alfaRad, sigms2, sigr1);
                    double p = pst + ptr;
                    p = Math.Round(p, 0);

                    //конец вытяжки с утонением
                    double msk = s[i] * Math.Sqrt(md[i]) / s[i - 1];
                    double sigms2k = sigmas2(sigmb, md[i], msk, psir);
                    sigms2k = Math.Round(sigms2k, 1);
                    double sigmr1k = 0;
                    double sigmzk = sigmaz2(msk, alfaRad, sigr1, sigms2k);
                    sigmzk = Math.Round(sigmzk, 1);
                    double pstk = pste(dsr[i], s[i], sigmzk);
                    double ptrk = Pi * ptre(msk, dsr[i], s[i], alfaRad, sigms2k, sigr1);
                    double pk = pstk + ptrk;
                    pk = Math.Round(pk, 0);
                    Console.WriteLine("p = " + p);
                    Console.WriteLine("sigmz = " + sigmz);
                    Console.WriteLine("sigms2 = " + sigms2);
                    Console.WriteLine("pk = " + pk);
                    Console.WriteLine("sigmzk = " + sigmzk);
                    Console.WriteLine("sigms2k = " + sigms2k);
                }
                else if (ms[i] == 1)
                {
                    double psisr = 1 - Math.Pow(md[i], 2);
                    double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
                    double sigmr = 1.1 * sigms * sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);
                    double p = Pi * dsr[i] * s[i] * sigmr;
                    sigms = Math.Round(sigms, 1);
                    sigmr = Math.Round(sigmr, 1);
                    p = Math.Round(p, 0);
                    Console.WriteLine("p = " + p);
                    Console.WriteLine("sigmr = " + sigmr);
                    Console.WriteLine("sigms = " + sigms);
                }
            }
        }

        public static double sigmar1(double alf, double md, double s, double dsr)
        {
            double mum = 0.05;
            return (1 + mum / Math.Tan(alf)) * Math.Log(1 / md) + 0.66 * Math.Sin(alf) * Math.Sqrt(s / dsr);
        }

        public static double sigmaz2(double ms, double alf, double sigmr1, double sigms2)
        {
            double mum = 0.05;
            double mup = 0.05;
            double ks = 1 / ms;
            double a = 1 - 0.5 * Math.Log(ks) - sigmr1;
            double b = 0.5 * (Math.Log(ks)) * (Math.Log(ks)) - sigmr1 * (ks - 1 - Math.Log(ks));
            return 1.15 * sigms2 * (Math.Log(ks) + mum * a * Math.Log(ks) / Math.Sin(alf) - mup * b / Math.Sin(alf) + sigmr1 + (1 - Math.Cos(alf)) / Math.Sin(alf));
        }

        public static double pste(double d, double s, double sigmz)
        {
            double Pi = Math.PI;
            return Pi * d * s * sigmz;
        }

        public static double ptre(double ms, double d, double s, double alf, double sigms2, double sigmr1)
        {
            double ks = 1 / ms;
            double a = Math.Log(ks) - sigmr1 * (ks - 1);
            return 1.15 * 0.05 * sigms2 * d * s * a / Math.Sin(alf);
        }

        public static double sigmas1(double sigmb, double md, double psir)
        {
            return sigmb * Math.Pow(((1 - md) / psir), (psir / (1 - psir)));
        }

        public static double D0KolpShtBezUt(double narDiam, double dopuNarZiam, double vnutrDiam, double dopuVnutrDiam, double vysota, double tolshchLenty, double verhDopTolst, double nizDopTolst, double radius)
        {
            double dn = narDiam + dopuNarZiam / 2;
            double dw = vnutrDiam + dopuVnutrDiam / 2;
            double h = vysota;
            double s = tolshchLenty + (verhDopTolst + nizDopTolst) / 2;
            double r = radius;
            double Pi = Math.PI;
            double vcil = Pi * (dn + dw) * s * (h - s - r) / 2;
            double vdno = Pi * s * (dw - 2 * r) * (dw - 2 * r) / 4;
            double vtor = Pi * Pi * (dw / 2 - r) * ((r + s) * (r + s) - r * r) / 2;
            double v = vcil + vdno + vtor;
            double d0KolpShtBezUt = Math.Sqrt(4 * v / Pi / s);
            return Math.Round(d0KolpShtBezUt, 2);
        }

        static void вырСверт2матр()
        {
            // Get the values from the user.
            Console.WriteLine("Enter the value of md1:");
            double md1 = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter the value of sd:");
            double sd = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter the value of psir:");
            double psir = Convert.ToDouble(Console.ReadLine());

            // Check if the shear stress is greater than (1 - md1) / 18.
            if (sd > (1 - md1) / 18)
            {
                Console.WriteLine("прижим не нужен на кон. и рад. матр");
            }
            else if (sd > (1 - md1) / 36)
            {
                Console.WriteLine("на конич. матр. прижим не нужен");
            }
            else
            {
                Console.WriteLine("нужен прижим");
            }

            // Get the type of matrix from the user.
            Console.WriteLine("Enter the type of matrix (матрица конич.-K, радиальная-R):");
            string tipMatriz = Console.ReadLine();

            // Check if the matrix type is conical.
            if (tipMatriz == "k")
            {
                if (sd > 0.012 && sd < 0.05)
                {
                    Console.WriteLine("Рекомендуется двухконусная матрица");
                }
                else if (sd > 0.05)
                {
                    Console.WriteLine("Рекомендуется одноконусная матрица с большим радиусом");
                }

                // Check if the shear stress is less than or equal to (1 - md1) / 36.
                if (sd <= (1 - md1) / 36)
                {
                    одноконМатрПрижим2матр();
                }
                else
                {
                    Console.WriteLine("матрица однокон.-O, двухкон.-D, однокон. с радиус-OR");
                    string tipKonMatriz = Console.ReadLine();
                    if (tipKonMatriz == "o")
                    {
                        одноконМатрБезПриж2матр();
                    }
                    else if (tipKonMatriz == "d")
                    {
                        двухконМатр2матр();
                    }
                    else if (tipKonMatriz == "or")
                    {
                        одноконМатрРадиус2матр();
                    }
                }
            }
            // Check if the matrix type is radial.
            else if (tipMatriz == "r")
            {
                if (sd > (1 - md1) / 18)
                {
                    радМатрБезПриж2матр();
                }
                else
                {
                    радМатрПрижим2матр();
                }
            }
        }

        public static void одноконМатрПрижим2матр()
        {
            // Get the input values.
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Get the user input for the angle of cone of the upper matrix.
            Console.WriteLine("Введи угол конусности верх. матрицы");
            double alf = double.Parse(Console.ReadLine());
            double al = Pi * alf / 180;

            double ms1 = double.Parse(Console.ReadLine());
            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the coefficient of friction.
            double mum = 0.05;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(Math.PI * alf / 180) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);
            double rw = 3 * s;
            double rws = rw + s / 2;

            // Calculate the force on the upper matrix
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double psi = 1 - dk / d0;
            double fik = Math.PI * (90 - alf) / 180;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double a = Math.Log(1 / md11) - psi + s / (2 * rws * Math.Pow(md11, 2));
            double b = 1 + mum * fik;
            double c = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * (a * b / (1 - 0.2 * mum * b * c / md1) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * (psisrn / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1sr == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s1 * sigmzk;
                //усилие при свертке с утонением в конечн. стадии;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Одноконус.с приж.через 2 матр.");
            double ustKonus = (Math.Sqrt((20 * sd) * (20 * sd) * (1 - Math.Sin(al)) + Math.Sin(al) * Math.Sin(al)) - 20 * sd) / Math.Sin(al);
            //проверка на склакообразованиие в конусе матрицы
            if (md1 < ustKonus)
            {
                Console.WriteLine("Однокон.с плоск.и КОНИЧ. приж. через 2 матр.");
            }

            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf);
        }

        public static void одноконМатрБезПриж2матр()
        {
            // Get the input values.
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Get the user input for the angle of cone of the upper matrix.
            Console.WriteLine("Введи угол конусности верх. матрицы");
            double alf = double.Parse(Console.ReadLine());
            double al = Pi * alf / 180;

            double ms1 = double.Parse(Console.ReadLine());
            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the coefficient of friction.
            double mum = 0.05;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = 0.9 * d0;
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);
            double rw = 0.5 * d0;
            double rws = rw + s / 2;
            double dkk = d1 * Math.Sqrt((1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(Math.PI * alf / 180) + 1);


            // Calculate the force on the upper matrix
            double md12 = d1 / dkk;
            double md11 = dkk / d0;
            double psi = 1 - dkk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);

            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(al)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s1 * sigmzk;
                //усилие при свертке с утонением в конечн. стадии;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Одноконус.матр.без приж.через 2 матр");
            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf);
        }

        public static void двухконМатр2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfanm = Math.Asin(q);
            double alfnm = alfanm * 180 / Pi;
            alfnm = Math.Round(alfnm, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfanm) + rn * (1 - Math.Sin(alfanm));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfanm)) * Math.Tan(alfanm);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //геометрия верхней матрицы
            double sd = double.Parse(Console.ReadLine());
            double werhugol;
            if (sd > 0.012 && sd < 0.018)
            {
                werhugol = 30;
                Console.WriteLine("Рекомендуемый угол верх. конуса = " + werhugol);
            }
            else if (sd > 0.0018 && sd < 0.05)
            {
                werhugol = 45;
                Console.WriteLine("Рекомендуемый угол верх. конуса = " + werhugol);
            }
            Console.WriteLine("Введи угол верх. конуса");
            double alfw = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи угол нижнего. конуса");
            double alfn = double.Parse(Console.ReadLine());
            double alw = Pi * alfw / 180;
            double alf = Pi * alfn / 180;

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfanm) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double dw = 0.9 * d0;
            dw = Math.Round(dw, 1);
            double hw = (dw - dk) / (2 * Math.Tan(alw));
            hw = Math.Round(hw, 1);
            double rw = 0.5 * d0;
            rw = Math.Round(rw, 1);
            double rs = (d0 - dm1) / 3;
            rs = Math.Round(rs, 1);
            double r = (d0 - d1) / 5;

            // Calculate the force on the upper matrix
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double fik = (alw + alf) / 2;
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (((1 + mum / Math.Tan(alw)) * (Math.Log(1 / md11) - psi) + s / (2 * r * Math.Sqrt(md11))) * (1 + mum * fik) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfanm / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfanm);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfnm);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Двухкон.матр2 матр");
            Console.WriteLine("Диам. верхнего конуса = " + dw);
            Console.WriteLine("Диаметр нижнего конуса = " + dk);
            Console.WriteLine("Высота верхнего конуса = " + hw);
            Console.WriteLine("Высота нижнего конуса = " + hk);
            Console.WriteLine("Угол верхнего конуса = " + alfw);
            Console.WriteLine("Угол верхнего конуса = " + alfn);
            Console.WriteLine("Радиус входн.верх конуса = " + rw);
            Console.WriteLine("Радиус входн.верх конуса = " + rs);
        }

        public static void одноконМатрРадиус2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            Console.WriteLine("Введи угол конуса верх. матрицы");
            double alfa = double.Parse(Console.ReadLine());
            double alf = Pi * alfa / 180;

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfa * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double md11 = dk / d0;
            double a = md11 + sd;
            double b = (1 - Math.Sin(alf)) * Math.Tan(alf);
            double hm = (a * (1 - b) + b - sd * ms1 - md1) * d0 / (2 * Math.Tan(alf));
            hm = Math.Round(hm, 1);
            double c = hm - hk;
            if (c < 0)
            {
                Console.WriteLine("высота матр. меньше высоты конуса");
                return;
            }
            double dkm = dm1 + 2 * hm * Math.Tan(alf);
            dkm = Math.Round(dkm, 1);
            double md12 = d1 / dk;
            double rw = (d0 - dk - s) / 2;
            rw = Math.Round(rw, 1);
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(alf)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Pi * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Однокон.матр.с больш.радиусом через 2матр");
            Console.WriteLine("Диам. конуса = " + dkm);
            Console.WriteLine("Высота конуса = " + hm);
            Console.WriteLine("Угол  конуса = " + alfa);
            Console.WriteLine("Радиус конуса = " + rw);
        }

        public static void радМатрПрижим2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());
            //нижняя матрица
            double ms1n = ms1 / ms1w;
            ms1n = Math.Round(ms1n, 3);

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;
            double dp = double.Parse(Console.ReadLine());
            d1 = dp + s1;
            //диаметр верхней матрицы
            double dm1 = dp + 2 * s1;
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;
            // угол кон. ниж. матр
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //верхней матрицы
            double rms = 0.69 * d1 * (0.16 * Math.Sqrt(18.3 / md1 * md1 + 21.2 + 10.2 * (rps / d1) * (rps / d1) - 41.3 * rps / d1) - 1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);
            double md11 = (d1 + 2 * rms) / d0;
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double fi = Math.PI / 2 - alr;
            double al = alr / 2;

            //расстояние между матрицами
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret;
            if (ms1w == 1)
            {
                tret = 0;
            }
            else
            {
                tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            }
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double a = Math.Log(1 / md1) + s / (4 * rm);
            double b = 1 + 1.5 * mum;
            double c = 0.2 * mum * b / md1;
            double d = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * a * b / (1 - c * d);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow(psisrn / psir, psir / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);

            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Радиал. матрица с прижимом 2 матр");
            if (sd < ((md11 - md1) / 20))
            {
                Console.WriteLine("Треб.плоск.и тор.прижим");
            }
            else
            {
                Console.WriteLine("Треб.плоск.и тор.прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }

        public static void радМатрБезПриж2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            //нижняя матрица
            double ms1n = ms1 / ms1w;
            ms1n = Math.Round(ms1n, 3);
            ms1 = ms1w;
            double s1 = s * ms1;
            double dp = double.Parse(Console.ReadLine());
            d1 = dp + s1;
            //диаметр верхней матрицы
            double dm1 = dp + 2 * s1;
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;
            // угол кон. ниж. матр
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //верхней матрицы
            double rms = d1 * (1 - md1) / (2 * md1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);
            double fig = 200 * (1 - md1 * (1 + 2.28 * rps / d1)) / (1 + md1 * (5 - 6 * md1));
            double fi = fig * Pi / 180;
            double md11 = md1 * (1 + 2 * rms * (1 - Math.Cos(fi)) / d1);
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double alf = alr / 2;

            //расстояние между матрицами
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret;
            if (ms1w == 1)
            {
                tret = 0;
            }
            else
            {
                tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            }
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum * fig) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow(psisrn / psir, psir / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);

            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Радиал. матрица без прижима 2матр");
            Console.WriteLine("Радиус матрицы = " + rm);
        }

        public static void koeffUton2матр()
        {
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the diameters and coefficients
            int[] dp = new int[n - 1];
            int[] dm = new int[n - 1];
            float[] msw = new float[n - 1];
            float[] msn = new float[n - 1];
            float[] ms = new float[n];

            // Read the data from the console
            for (int i = 0; i < n - 1; i++)
            {
                dp[i] = int.Parse(Console.ReadLine());
                dm[i] = int.Parse(Console.ReadLine());
            }

            // Calculate the total coefficient of attenuation for each operation
            for (int j = 0; j < n; j++)
            {
                ms[j] = msw[j] * msn[j];
            }
        }
        public double AngleOfRadialMatrix(double previousDiameter, double diameter, double radius)
        {
            double d1 = previousDiameter;
            double d2 = diameter;
            double r = radius;
            double a = r - (d1 - d2) / 2;
            double b = Math.Sqrt(r * r - a * a);
            double c = (d1 - d2) / 2;
            double alfa = Math.Atan(c / b);
            alfa *= 180 / 3.1416;
            return alfa;
        }

        public double ThicknessOfConeWall(double drawRatio, double coneAngle, double frictionCoefficient, double initialThickness)
        {
            double md = drawRatio;
            double alfa = coneAngle;
            double Pi = 3.1416;
            alfa *= Pi / 180;
            double mu = frictionCoefficient;
            double s0 = initialThickness;
            double a = 1 + mu / Math.Tan(alfa);
            double b = 1 / md - 1;
            double c = 0.5 - 0.75 * a * b / (2 - 0.5 * a * b);
            double s = s0 * Math.Pow(1 / md, c);
            return s;
        }

        public void ВырВытяжкаСдвигом()
        {
            double dn = double.Parse(Console.ReadLine());
            double deltan = double.Parse(Console.ReadLine());
            double dw = double.Parse(Console.ReadLine());
            double deltaw = double.Parse(Console.ReadLine());
            double d1 = (dn + deltan / 2 + dw + deltaw) / 2;
            double s1 = ((dn + deltan / 2) - (dw + deltaw));
            double s0 = double.Parse(Console.ReadLine());
            double sd1 = s1 / d1;
            Console.WriteLine("Предельный коэф. утонения для стали и алюмю=0,4 для латуни=0,35");
            Console.WriteLine("Введи угол конусности матрицы");
            double alf = double.Parse(Console.ReadLine());
            Console.WriteLine(sd1);
            double kt1 = s0 / s1;
            double al = 3.1416 * alf / 180;
            double md1 = 1 / (1 + sd1 * (2 - 1 / kt1));
            double m = 1 + (1 - md1 * md1) * kt1 * kt1 * Math.Tan(al) / (2 * sd1 * md1 * md1);
            double ms1 = 1 / Math.Sqrt(m);
            Console.WriteLine(ms1);
            Console.WriteLine("штамп. сдвигом");
        }
    }
}

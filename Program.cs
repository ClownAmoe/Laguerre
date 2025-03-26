using System;

public class LaguerreFunc
{
    private int _n;
    private float _t;
    private float _sigma;
    private float _beta;

    public LaguerreFunc(int n, double t, double sigma = 4, double beta = 2)
    {
        N = n;
        T = (float)t;
        Sigma = (float)sigma;
        Beta = (float)beta;
    }

    public int N
    {
        get { return _n; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("N cannot be negative.");
            }
            _n = value;
        }
    }

    public float T
    {
        get { return _t; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("T cannot be negative.");
            }
            _t = value;
        }
    }

    public float Sigma
    {
        get { return _sigma; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("Sigma cannot be negative.");
            }
            _sigma = value;
        }
    }

    public float Beta
    {
        get { return _beta; }
        set
        {
            if (value < 0)
            {
                throw new ArgumentException("Beta cannot be negative.");
            }
            _beta = value;
        }
    }

    public double CalcFunc()
    {
        double l0 = Math.Sqrt(Sigma) * Math.Exp(-Beta / 2 * T);
        double l1 = Math.Sqrt(Sigma) * (1 - Sigma * T) * Math.Exp(-Beta / 2 * T);

        if (_n == 0)
        {
            return l0;
        }
        else if (_n == 1)
        {
            return l1;
        }

        double ln2 = l0;
        double ln1 = l1;
        double ln = 0;

        for (int i = 2; i <= _n; ++i)
        {
            ln = ((2 * i - 1 - Sigma * T) / i) * ln1 - ((i - 1) / i) * ln2;
            ln2 = ln1;
            ln1 = ln;
        }

        return ln;
    }

    public (double[], double[]) TabulateLaguerre(double T, int numPoints = 100)
    {
        double[] tValues = new double[numPoints];
        double[] lValues = new double[numPoints];

        for (int i = 0; i < numPoints; i++)
        {
            tValues[i] = i * (T / (numPoints - 1));
            this.T = (float)tValues[i];
            lValues[i] = CalcFunc();
        }

        return (tValues, lValues);
    }

    public double LeftRect(double alpha, double T, double epsilon = 1e-3)
    {
        int numPoints = (int)Math.Ceiling(T / epsilon);
        double dt = T / numPoints;
        double integralSum = 0;

        for (int i = 0; i < numPoints; i++)
        {
            double t = i * dt;
            this.T = (float)t;
            integralSum += CalcFunc() * Math.Exp(-alpha * t) * dt;
        }

        return integralSum;
    }

    public double MidRect(double alpha, double T, double epsilon = 1e-3)
    {
        int numPoints = (int)Math.Ceiling(T / epsilon);
        double dt = T / numPoints;
        double integralSum = 0;

        for (int i = 0; i < numPoints; i++)
        {
            double t = (i + 0.5) * dt;
            this.T = (float)t;
            integralSum += CalcFunc() * Math.Exp(-alpha * t) * dt;
        }

        return integralSum;
    }

    public double InverseLaguerre(double[] h)
    {
        int temp = _n;
        double sum = 0;

        for (int i = 0; i < h.Length; ++i)
        {
            _n = i;
            sum += h[i] * CalcFunc();
        }

        _n = temp;
        return sum;
    }
    public double FindT(int N = 20, double epsilon = 1e-3)
    {
        double T = 0.0;
        bool found = false;

        while (!found)
        {
            found = true;

            for (int n = 0; n <= N; n++)
            {
                this.N = n;
                this.T = (float)T;

                if (Math.Abs(CalcFunc()) >= epsilon)
                {
                    found = false;
                    break;
                }
            }

            if (!found)
            {
                T += 0.1;
            }
        }

        return T;
    }
}


class ProgramLaguerre
{
    static void Main()
    {
        double[] h = { 1, 2, 3, 4, 5 };
        var lag = new LaguerreFunc(0, 1);
        Console.WriteLine("Laguerre value: " + lag.CalcFunc());
        Console.WriteLine("Inverse: " + lag.InverseLaguerre(h));

        Console.WriteLine("Left rect: " + lag.LeftRect(1, 5));
        Console.WriteLine("Mid rect: " + lag.MidRect(1, 5));

        double T = lag.FindT();
        Console.WriteLine($"T = {T}");


        var (tValues, lValues) = lag.TabulateLaguerre(T);
        Console.WriteLine("t\t| L(n, t)");
        Console.WriteLine("----------------");

        for (int i = 0; i < tValues.Length; i++)
        {
            Console.WriteLine($"{tValues[i]:F3}\t| {lValues[i]:F6}");
        }
    }
}
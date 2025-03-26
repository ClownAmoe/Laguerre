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


    public double InverseLaguerre(double[] h)
    {
        int temp = _n;
        double sum = 0;

        for (int i = 0; i < _n; ++i)
        {
            _n = i;
            sum += h[i] * CalcFunc();
        }

        _n = temp;
        return sum;
    }
}

class ProgramLaguerre
{
    static void Main()
    {
        double[] h = { 1, 2, 3, 4, 5 };
        var lag = new LaguerreFunc(0, 1);
        Console.WriteLine(lag.CalcFunc());
        Console.WriteLine(lag.InverseLaguerre(h));
    }
}
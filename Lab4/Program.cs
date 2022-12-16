using ScottPlot;
internal class Program
{
    struct Spectrum
    {
        public double[] amplitude;
        public double[] phase;

        public Spectrum(double[] amplitude, double[] phase)
        {
            this.amplitude = amplitude;
            this.phase = phase;
        }
    }

    private static void Main(string[] args)
    {
        int width = 1920;
        int height = 1080;
        var N = 1024;
        var count = new double[N];
        for (int i = 0; i < N; i++)
        {
            count[i] = i;
        }
        var signal = RandomSignal(N);
        var rollingAvgValues = RollingAverage(signal, 3);
        var parabolaValues = Parabola(signal);
        var rollingMedianValues = RollingMedian(signal, 5);

        //Task 3
        var plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        plot.AddScatter(count, signal, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, rollingAvgValues, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Smoothed by moving average");
        plot.SaveFig("Task3-1.png");
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        plot.AddScatter(count, signal, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, parabolaValues, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Quarter-smoothed parabola (7 points)");
        plot.SaveFig("Task3-2.png");
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        plot.AddScatter(count, signal, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, rollingMedianValues, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Smoothed by median filtering");
        plot.SaveFig("Task3-3.png");

        //Task 5
        var signalSpectrum = GetSpectrum(signal);
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        var sigPhase = plot.AddSignal(signalSpectrum.phase, label: "PS");
        sigPhase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        var sigAmp = plot.AddSignal(signalSpectrum.amplitude, label: "AS");
        var yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-1.png");

        var rollingAvgAASpectrum = GetSpectrum(rollingAvgValues);
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        sigPhase = plot.AddSignal(rollingAvgAASpectrum.phase, label: "PS");
        sigPhase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(rollingAvgAASpectrum.amplitude, label: "AS");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-2.png");

        var parabolaAASpectrum = GetSpectrum(parabolaValues);
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        sigPhase = plot.AddSignal(parabolaAASpectrum.phase, label: "PS");
        sigPhase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(parabolaAASpectrum.amplitude, label: "AS");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-3.png");

        var rollingMedianAASpectrum = GetSpectrum(rollingMedianValues);
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        sigPhase = plot.AddSignal(rollingMedianAASpectrum.phase, label: "PS");
        sigPhase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(rollingMedianAASpectrum.amplitude, label: "AS");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-4.png");
    }

    static Spectrum GetSpectrum(double[] signal)
    {
        var N = signal.Length;
        Spectrum result = new Spectrum(new double[N], new double[N]);
        Func<double[], Func<double, double>, int, double> Dftmplitude =
            (signal, MathFunc, j) =>
            {
                double resultAmp = 0;
                for (int i = 0; i < signal.Length; i++)
                {
                    resultAmp += signal[i] * MathFunc(2 * Math.PI * j * i / signal.Length);
                }
                return 2 * resultAmp / signal.Length;
            };

        for (int j = 0; j < N; j++)
        {
            double cos = Dftmplitude(signal, Math.Cos, j);
            double sin = Dftmplitude(signal, Math.Sin, j);
            result.amplitude[j] = Math.Sqrt(Math.Pow(cos, 2) + Math.Pow(sin, 2));
            result.phase[j] = Math.Atan2(sin, cos);
        }
        return result;
    }

    static double[] RollingAverage(double[] signal, int windowSize)
    {
        var N = signal.Length;
        double[] result = new double[N];
        int M = (windowSize - 1) / 2;
        for (int i = 0; i < N; i++)
        {
            int low = (i - M) < 0 ? 0 : (i - M);
            int high;
            if (low == 0)
            {
                high = low + windowSize - 1;
            }
            else if ((i + M) >= N)
            {
                high = N - 1;
                low = high - windowSize + 1;
            }
            else
            {
                high = i + M;
            }
            for (int j = low; j <= high; j++)
            {
                result[i] += signal[j];
            }
            result[i] /= windowSize;
        }
        return result;
    }

    static double[] Parabola(double[] signal)
    {
        var N = signal.Length;
        double[] result = new double[N];
        var windowSize = 7;
        var divisor = 231;
        double[] multipliers = new double[] { 5, -30, 75, 131, 75, -30, 5 };
        int M = (windowSize - 1) / 2;
        for (int i = 0; i < N; i++)
        {
            int low = (i - M) < 0 ? 0 : (i - M);
            int high;
            if (low == 0)
            {
                high = low + windowSize - 1;
            }
            else if ((i + M) >= N)
            {
                high = N - 1;
                low = high - windowSize + 1;
            }
            else
            {
                high = i + M;
            }

            for (int j = low; j <= high; j++)
            {
                result[i] += multipliers[j - low] * signal[j];
            }
            result[i] /= divisor;
        }
        return result;
    }

    static double[] RollingMedian(double[] signal, int windowSize)
    {
        var N = signal.Length;
        double[] result = new double[N];
        double[] tmp = new double[windowSize];
        int K = 1;
        for (int i = 0; i < N; i++)
        {
            int low = (i < windowSize) ? 0 : i;
            if ((i + windowSize) >= N)
            {
                low = N - 1 - windowSize;
            }
            Array.Copy(signal, low, tmp, 0, windowSize);
            Array.Sort(tmp);
            for (int j = K; j < tmp.Length - K; j++)
            {
                result[i] += tmp[j];
            }
            result[i] /= (windowSize - 2 * K);
        }
        return result;
    }

    static double[] RandomSignal(int n)
    {
        var result = new double[n];
        Random rand = new Random(999999);
        Func<Random, int, int, double> noisy = (rand, i, n) =>
        {
            double value = 0;
            for (int j = 50; j <= 70; j++)
            {
                value += 2 * Math.Sin(2 * Math.PI * i * j / n) * (rand.NextDouble() < 0.5 ? -1 : 1);
            };
            return value;
        };

        for (int i = 0; i < n; i++)
        {
            result[i] = 4 * Math.Sin(2 * Math.PI * i / n) + noisy(rand, i, n);
        }
        return result;
    }
}
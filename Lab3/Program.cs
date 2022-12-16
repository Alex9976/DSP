using ScottPlot;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics;

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

    static double[] SinTable;
    static Func<int, bool> HighPassFilter = (i) => i > 20;
    static Func<int, bool> LowPassFilter = (i) => i < 30;
    static Func<int, bool> BandPassFilter = (i) => (i < 10 || i > 15);
    static Func<int, int, double> TableSin = (i, j) => SinTable[(i * j) % SinTable.Length];
    static Func<int, int, double> TableCos = (i, j) => SinTable[(i * j + SinTable.Length / 4) % SinTable.Length];

    private static void Main(string[] args)
    {
        int width = 1920;
        int height = 1080;
        var N = 1024;
        var count = new double[N];
        var amplitude = 10;
        double[] SignalValues = HarmonicSignal(amplitude, 0, 1, N);
        SinTable = InitSinTable(N, amplitude);
        var spectrum = TableSpectrum(SignalValues, amplitude);
        double[] restoredSignal = RestoreSignalFromSpectrum(spectrum);

        //Task 2ab
        var plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        for (int i = 0; i < N; i++)
        {
            count[i] = i;
        }
        plot.AddScatter(count, SignalValues, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, restoredSignal, markerShape: MarkerShape.none, lineStyle: LineStyle.Dot, label: "Restored");
        plot.SaveFig("Task2ab-1.png");

        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        var phase = plot.AddSignal(spectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        var sigAmp = plot.AddSignal(spectrum.amplitude, label: "Amplitude spectrum");
        var yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task2ab-2.png");

        //Task 3abc
        var amplitudes = new double[] { 1, 3, 5, 8, 10, 12, 16 };
        var phases = new double[] { Math.PI / 6, Math.PI / 4, Math.PI / 3, Math.PI / 2, 3 * Math.PI / 4, Math.PI };
        var polyharmonicSignal = PolyharmonicSignal(amplitudes, phases, N);
        var signalSpectrum = GetSpectrum(polyharmonicSignal);

        double[] restoredSignalValues = RestoreSignalFromSpectrum(spectrum);
        Spectrum signalSpectrumWOPhase = new Spectrum(spectrum.amplitude, new double[N]);
        double[] restoredSignalWOPhase = RestoreSignalFromSpectrum(signalSpectrumWOPhase);

        for (int i = 0; i < restoredSignal.Length; i++)
        {
            restoredSignal[i] += amplitudes[0] / 2;
            restoredSignalWOPhase[i] += amplitudes[0] / 2;
        }

        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        for (int i = 0; i < N; i++)
        {
            count[i] = i;
        }
        plot.AddScatter(count, polyharmonicSignal, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, restoredSignal, markerShape: MarkerShape.none, lineStyle: LineStyle.Dot, label: "Restored");
        plot.AddScatter(count, restoredSignalWOPhase, markerShape: MarkerShape.none, lineStyle: LineStyle.DashDot, label: "Restored without phase");
        plot.SaveFig("Task3abc-1.png");

        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        phase = plot.AddSignal(spectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(spectrum.amplitude, label: "Amplitude spectrum");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task3abc-2.png");

        // Task 4
        var fftSpectrum = Fft(polyharmonicSignal);
        double[] restoredFftSignal = RestoreSignalFromSpectrum(fftSpectrum);
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        plot.AddScatter(count, polyharmonicSignal, markerShape: MarkerShape.none, lineStyle: LineStyle.Solid, label: "Source");
        plot.AddScatter(count, restoredFftSignal, markerShape: MarkerShape.none, lineStyle: LineStyle.Dot, label: "Restored");
        plot.SaveFig("Task4-1.png");

        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        phase = plot.AddSignal(fftSpectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(fftSpectrum.amplitude, label: "Amplitude spectrum");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task4-2.png");

        var hpFilteredSpectrum = FilterSpectrum(spectrum, HighPassFilter);
        var lpFilteredSpectrum = FilterSpectrum(spectrum, LowPassFilter);
        var bpFilteredSpectrum = FilterSpectrum(spectrum, BandPassFilter);

        var hpSignal = RestoreSignalFromSpectrum(hpFilteredSpectrum);
        var lpSignal = RestoreSignalFromSpectrum(lpFilteredSpectrum);
        var bpSignal = RestoreSignalFromSpectrum(bpFilteredSpectrum);

        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        plot.AddScatter(count, restoredSignal, lineStyle: LineStyle.Solid, label: "Restored");
        plot.AddScatter(count, hpSignal, lineStyle: LineStyle.Dot, label: "Signal(HP)");
        plot.AddScatter(count, lpSignal, lineStyle: LineStyle.DashDot, label: "Signal(LP)");
        plot.AddScatter(count, bpSignal, lineStyle: LineStyle.DashDot, label: "Signal(BP)");
        plot.SaveFig("Task5.png");

        //HP
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        phase = plot.AddSignal(hpFilteredSpectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(hpFilteredSpectrum.amplitude, label: "Amplitude spectrum");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-HP.png");

        //LP
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        phase = plot.AddSignal(lpFilteredSpectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(lpFilteredSpectrum.amplitude, label: "Amplitude spectrum");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-LP.png");

        //BP
        plot = new Plot(width, height);
        plot.Legend(true, Alignment.UpperCenter);
        phase = plot.AddSignal(bpFilteredSpectrum.phase, label: "Phase spectrum");
        phase.YAxisIndex = 0;
        plot.YAxis.Label("PS");

        sigAmp = plot.AddSignal(bpFilteredSpectrum.amplitude, label: "Amplitude spectrum");
        yAxis3 = plot.AddAxis(ScottPlot.Renderable.Edge.Left, axisIndex: 2);
        sigAmp.YAxisIndex = 2;
        yAxis3.Label("AS");
        plot.SaveFig("Task5-BP.png");
    }

    //Task 4 FFT
    static Spectrum Fft(double[] signal)
    {
        Complex32[] data = new Complex32[signal.Length];
        for (int i = 0; i < signal.Length; i++)
        {
            data[i] = (float)signal[i];
        }

        Fourier.Forward(data, FourierOptions.NumericalRecipes);
        Spectrum result = new Spectrum(new double[signal.Length], new double[signal.Length]);
        for (int i = 0; i < signal.Length; i++)
        {
            data[i] = new Complex32(data[i].Real * 2 / signal.Length, data[i].Imaginary * 2 / signal.Length);
            result.amplitude[i] = Math.Sqrt(Math.Pow(data[i].Real, 2) + Math.Pow(data[i].Imaginary, 2));
            result.phase[i] = Math.Atan2((data[i].Imaginary), (data[i].Real));
        }
        return result;
    }

    //Task 5
    static Spectrum FilterSpectrum(Spectrum spectrum, Func<int, bool> filterFreq)
    {
        var N = spectrum.amplitude.Length;
        Spectrum result = new Spectrum(new double[N], new double[N]);
        Array.Copy(spectrum.amplitude, result.amplitude, N);
        Array.Copy(spectrum.phase, result.phase, N);
        var half = result.amplitude.Length / 2;
        for (int i = 0; i < half; i++)
        {
            if (filterFreq(i))
            {
                result.amplitude[i] = result.amplitude[N - 1 - i] = result.phase[i] = result.phase[N - 1 - i] = 0;
            }
        }
        return result;
    }

    static double[] InitSinTable(int N, int amplitude)
    {
        var result = new double[N];
        for (int i = 0; i < result.Length; i++)
        {
            result[i] = amplitude * Math.Sin(2 * Math.PI * i / N);
        }
        return result;
    }

    static Spectrum TableSpectrum(double[] signal, int A)
    {
        var N = signal.Length;
        Spectrum result = new Spectrum(new double[N], new double[N]);
        for (int j = 0; j < N; j++)
        {
            double cos = DftAmplitude(signal, TableCos, A, j);
            double sin = DftAmplitude(signal, TableSin, A, j);
            result.amplitude[j] = Math.Sqrt(Math.Pow(cos, 2) + Math.Pow(sin, 2));
            result.phase[j] = Math.Atan2(sin, cos);
        }
        return result;
    }

    static Spectrum GetSpectrum(double[] values)
    {
        var N = values.Length;
        Spectrum result = new Spectrum(new double[N], new double[N]);
        Func<double[], Func<double, double>, int, double> Dftmplitude =
            (signal, MathFunc, j) =>
            {
                double amplitude = 0.0;
                for (int i = 0; i < signal.Length; i++)
                {
                    amplitude += signal[i] * MathFunc(2 * Math.PI * j * i / signal.Length);
                }
                return 2 * amplitude / signal.Length;
            };

        for (int j = 0; j < N; j++)
        {
            double cos = Dftmplitude(values, Math.Cos, j);
            double sin = Dftmplitude(values, Math.Sin, j);
            result.amplitude[j] = Math.Sqrt(Math.Pow(cos, 2) + Math.Pow(sin, 2));
            result.phase[j] = Math.Atan2(sin, cos);
        }
        return result;
    }

    static double HarmonicCos(double amplitude, double phase, double freq, int n, int i)
    {
        return amplitude * Math.Cos((2 * Math.PI * freq * i) / n + phase);
    }

    static double RestoredSignal(double[] amplitude, double[] phase, int i)
    {
        double result = 0;
        for (int j = 0; j < amplitude.Length / 2; j++)
        {
            result += amplitude[j] * Math.Cos((2 * Math.PI * j * i / amplitude.Length) - phase[j]);
        }
        return result;
    }

    static double[] PolyharmonicSignal(double[] amplitude, double[] phases, int n)
    {
        var result = new double[n];
        Random rand = new Random();
        Func<double[], double[], Random, int, int, double> value = (amplitudes, phases, rand, n, i) =>
        {
            double result = 0;
            rand = new Random(999999);
            for (int j = 1; j <= 30; j++)
            {
                result += amplitudes[rand.Next(amplitudes.Length - 1)] * Math.Cos(2 * Math.PI * j * i / n - phases[rand.Next(phases.Length - 1)]);
            }
            return result;
        };

        for (int i = 0; i < n; i++)
        {
            result[i] = value(amplitude, phases, rand, n, i);
        }
        return result;
    }

    static double[] HarmonicSignal(double amplitude, double phase, double freq, int n)
    {
        var result = new double[n];
        Parallel.For(0, n, i =>
        {
            result[i] = HarmonicCos(amplitude, phase, freq, n, i);
        });
        return result;
    }

    static double DftAmplitude(double[] signal, Func<int, int, double> table, int amplitude, int j)
    {
        double resultAmp = 0.0;
        for (int i = 0; i < signal.Length; i++)
        {
            resultAmp += signal[i] * table(i, j);
        }
        return 2 * resultAmp / (amplitude * signal.Length);
    }

    static double[] RestoreSignalFromSpectrum(Spectrum spectrum)
    {
        var result = new double[spectrum.amplitude.Length];
        for (int i = 0; i < result.Length; i++)
        {
            result[i] = RestoredSignal(spectrum.amplitude, spectrum.phase, i);
        }
        return result;
    }
}
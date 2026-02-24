

using ConsoleGPS;

class Program
{
    static void Main()
    {
        Console.WriteLine("Преобразование координат демонстрационный вариант");

        // Широта Долгота исходные
        string sB = "53°12'05.41";
        string sL = "50°06'48.03";

        double B = ParserGMS.ParseToDecimal(sB) * Math.PI / 180.0;
        double L = ParserGMS.ParseToDecimal(sL) * Math.PI / 180.0;

        double[] o1 = TranslationWGS84.ConvertToGaussKrueger(B, L);                 // Плоские координаты для исходных данных
        double[] testBl = TranslationWGS84.ConvertFromGaussKrueger(o1[0], o1[1]);   // Широта Долгота полученные

        Console.WriteLine("Исходные широта долгота:");
        Console.WriteLine(sB + " " + sL + "\n");

        Console.WriteLine("Пересчитанные широта долгота:");
        Console.WriteLine(ParserGMS.ParseToDMS(testBl[0] * 180 / Math.PI) + " " + ParserGMS.ParseToDMS(testBl[1] * 180 / Math.PI) + "\n");


        double[] o2 = TranslationWGS84.ConvertToGaussKrueger(testBl[0], testBl[1]);
        Console.WriteLine();
        Console.WriteLine("Сравним плоские координаты");
        Console.WriteLine(o1[0] + " " + o1[1]);
        Console.WriteLine(o2[0] + " " + o2[1]);

        Console.WriteLine("Разница между исходной точкой и повторным её пересчётом: " +
                Math.Sqrt(Math.Pow(o1[0] - o2[0], 2) + Math.Pow(o1[1] - o2[1], 2)));
    }
}


/*
        double[] d_PZ90 = TranslationWGS84.ConvertCartesianToGeodetic(
            TranslationWGS84.ConvertWGS84ToPZ90(TranslationWGS84.ConvertGeodeticToCartesian(0.523599, 0.523599, EncodingTypeGPS.WGS_84))
            , EncodingTypeGPS.PZ_90);
        Console.WriteLine("B = " + d_PZ90[0]);
        Console.WriteLine("L = " + d_PZ90[1]);
        Console.WriteLine("" + Math.Pow(10, -8));

        Console.WriteLine();
        double[] test = TranslationWGS84.ConvertCartesianToGeodetic(TranslationWGS84.ConvertGeodeticToCartesian(0.035, 0.066, EncodingTypeGPS.WGS_84), EncodingTypeGPS.WGS_84);

        Console.WriteLine(test[0] + " " + test[1] + " " + test[2]);


        Console.WriteLine();
        double[] o1 = TranslationWGS84.ConvertToGaussKrueger(
            Math.PI * (53.0 + 12.0 / 60.0 + 59.53 / 3600.0) / 180.0,
            Math.PI * (50.0 + 7.0 / 60.0 + 54.32 / 3600.0) / 180.0
            );

        double[] o2 = TranslationWGS84.ConvertToGaussKrueger(
            Math.PI * (53.0 + 12.0 / 60.0 + 58.62 / 3600.0) / 180.0,
            Math.PI * (50.0 + 7.0 / 60.0 + 52.49 / 3600.0) / 180.0
            );

        Console.WriteLine("" + Math.Sqrt(Math.Pow(o1[0] - o2[0], 2) + Math.Pow(o1[1] - o2[1], 2)));
        */
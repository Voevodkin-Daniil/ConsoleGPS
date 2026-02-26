

using ConsoleGPS;

class Program
{
    static void Main()
    {
        Console.WriteLine("Преобразование координат демонстрационный вариант");

        // Широта Долгота исходные
        Console.WriteLine("Преобразование координат демонстрационный вариант");


        // Сравнение координат на плоскости для БПЛА
        double BPLA_B1 = 53.2017412 * Math.PI / 180.0;
        double BPLA_L1 = 50.1117721 * Math.PI / 180.0;

        double BPLA_B2 = 53.2017412663672 * Math.PI / 180.0;
        double BPLA_L2 = 50.11177213010242 * Math.PI / 180.0;

        var BPLA_coords1 = TranslationWGS84.ConvertToGaussKrueger(BPLA_B1, BPLA_L1);
        var BPLA_coords2 = TranslationWGS84.ConvertToGaussKrueger(BPLA_B2, BPLA_L2);

        Console.WriteLine("Разница расстояний между истинным и найденным для БПЛА");
        Console.WriteLine(Math.Sqrt(Math.Pow(BPLA_coords1[0] - BPLA_coords2[0], 2) +
                                    Math.Pow(BPLA_coords1[1] - BPLA_coords2[1], 2)));

        Console.WriteLine();

        // Сравнение координат на плоскости для Цели
        double target_B1 = 53.2007990 * Math.PI / 180.0;
        double target_L1 = 50.1142422 * Math.PI / 180.0;

        double target_B2 = 53.20079887330737 * Math.PI / 180.0;
        double target_L2 = 50.11424239584601 * Math.PI / 180.0;

        var target_coords1 = TranslationWGS84.ConvertToGaussKrueger(target_B1, target_L1);
        var target_coords2 = TranslationWGS84.ConvertToGaussKrueger(target_B2, target_L2);

        Console.WriteLine("Разница расстояний между истинным и найденным для ЦЕЛИ");
        Console.WriteLine(Math.Sqrt(Math.Pow(target_coords1[0] - target_coords2[0], 2) +
                                    Math.Pow(target_coords1[1] - target_coords2[1], 2)));
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
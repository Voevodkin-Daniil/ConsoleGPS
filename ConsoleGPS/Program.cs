

using ConsoleGPS;

class Program
{
    static void Main()
    {
        Console.WriteLine("Преобразование координат демонстрационный вариант");


        // Сравнение координат на плоскости для БПЛА
        double BPLA_B1 = 53.2017412 * Math.PI / 180.0;
        double BPLA_L1 = 50.1117721 * Math.PI / 180.0;

        double BPLA_B2 = 53.2017412663672 * Math.PI / 180.0;
        double BPLA_L2 = 50.11177213010242 * Math.PI / 180.0;

        var BPLA_coords1 = TranslationWGS84.ConvertToGaussKrueger(BPLA_B1, BPLA_L1);
        var BPLA_coords2 = TranslationWGS84.ConvertToGaussKrueger(BPLA_B2, BPLA_L2);

        Console.WriteLine("Разница расстояний между истинным и найденным для БПЛА");
        Console.WriteLine(Math.Sqrt(Math.Pow(BPLA_coords1.x - BPLA_coords2.x, 2) +
                                    Math.Pow(BPLA_coords1.y - BPLA_coords2.y, 2)));

        Console.WriteLine();

        // Сравнение координат на плоскости для Цели
        double target_B1 = 53.2007990 * Math.PI / 180.0;
        double target_L1 = 50.1142422 * Math.PI / 180.0;

        double target_B2 = 53.20079887330737 * Math.PI / 180.0;
        double target_L2 = 50.11424239584601 * Math.PI / 180.0;

        var target_coords1 = TranslationWGS84.ConvertToGaussKrueger(target_B1, target_L1);
        var target_coords2 = TranslationWGS84.ConvertToGaussKrueger(target_B2, target_L2);

        Console.WriteLine("Разница расстояний между истинным и найденным для ЦЕЛИ");
        Console.WriteLine(Math.Sqrt(Math.Pow(target_coords1.x - target_coords2.x, 2) +
                                    Math.Pow(target_coords1.y - target_coords2.y, 2)));

    }
}

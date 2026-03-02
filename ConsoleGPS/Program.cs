

using ConsoleGPS;

class Program
{
    static void Main()
    {
        Console.WriteLine("Преобразование координат демонстрационный вариант");

        (double B, double L, double H) OO1 = (53.2015043 * Math.PI / 180.0, 50.1133423 * Math.PI / 180.0, 0);

        (double B, double L, double H) OO2 = (53.2000976 * Math.PI / 180.0, 50.1125830 * Math.PI / 180.0, 0);

        // Сравнение насколько будет разница при двойном определении
        test_1(OO1);

        // Проерка что выдаст растояние
        test_2(OO1, OO2);
        test_3(OO2, OO1);

        // Сравнение ответов при преобразовании
        test_4();

    }

    private static bool test_1((double B, double L, double H) Point)
    {
        var temp = TranslationWGS84.ConvertToGaussKrueger(Point.B, Point.L);
        var point_1 = TranslationWGS84.ConvertFromGaussKrueger(temp.x, temp.y);

        var temp2 = TranslationWGS84.ConvertToGaussKrueger(point_1.B, point_1.L);

        Console.WriteLine(Math.Sqrt(Math.Pow(temp.x - temp2.x,2) + Math.Pow(temp.y - temp2.y, 2)));
        Console.WriteLine();

        return true;
    }

    private static bool test_2((double B, double L, double H) A, (double B, double L, double H) B)
    {
        var temp_A = TranslationWGS84.ConvertToGaussKrueger(A.B, A.L);
        var temp_B = TranslationWGS84.ConvertToGaussKrueger(B.B, B.L);

        Console.WriteLine(Math.Sqrt(Math.Pow(temp_A.x - temp_B.x, 2) + Math.Pow(temp_A.y - temp_B.y, 2)));
        Console.WriteLine();

        return true;
    }

    private static bool test_3((double B, double L, double H) T, (double B, double L, double H) OO)
    {
        var temp_T = TranslationWGS84.ConvertGeodeticToLocalTangent(T, OO);

        
        Console.WriteLine(Math.Sqrt(Math.Pow(temp_T.east, 2) + Math.Pow(temp_T.north, 2)));

        var T_new = TranslationWGS84.ConvertLocalTangentToGeodetic(temp_T, OO);
        Console.WriteLine($"{(T.B / Math.PI * 180):f7} \t {(T.L / Math.PI * 180):f7}");
        Console.WriteLine($"{(T_new.B / Math.PI * 180):f7} \t {(T_new.L / Math.PI * 180):f7}");
        Console.WriteLine();

        return true;
    }

    private static bool test_4()
    {
        Console.WriteLine("Тест 4");
        // Начальные значения
        double OO1_B = 53.2015043 / 180 * Math.PI;
        double OO1_L = 50.1133423 / 180 * Math.PI;
        double OO1_H = -10;
        double OO2_B = 53.2000976 / 180 * Math.PI;
        double OO2_L = 50.1125830 / 180 * Math.PI;
        double OO2_H = 13;
        double L1 = 210.01;
        double L2 = 247.08;
        double L3 = 246.46;

        double a = -18.32;
        double b = 41.07;
        double aa = -58.99;
        double bb = -39.45;
        double cc = -37.48;


        // Проверка на проекции Гаусса-Крюгера
        var temp_OO1_1 = TranslationWGS84.ConvertToGaussKrueger(OO1_B, OO1_L);
        var temp_OO2_1 = TranslationWGS84.ConvertToGaussKrueger(OO2_B, OO2_L);

        (var B1, var T1) = TargetLocator.GetTargetCoordinates(
                            temp_OO1_1.x, temp_OO1_1.y, OO1_H,
                            temp_OO2_1.x, temp_OO2_1.y, OO2_H,
                            L1, L2, L3, a, b, aa, bb, cc);

        var gB1 = TranslationWGS84.ConvertFromGaussKrueger(B1.X, B1.Y);
        var gT1 = TranslationWGS84.ConvertFromGaussKrueger(T1.X, T1.Y);

        Console.WriteLine($"B = {(gB1.B / Math.PI * 180):f7}, {(gB1.L / Math.PI * 180):f7}, {B1.Z:f2}");
        Console.WriteLine($"T = {(gT1.B / Math.PI * 180):f7}, {(gT1.L / Math.PI * 180):f7}, {T1.Z:f2}");
        Console.WriteLine();

        // Проверка на Координатах локальной косательной плоскости

        var temp_OO1_2 = TranslationWGS84.ConvertGeodeticToLocalTangent(
            (OO1_B, OO1_L, OO1_H),
            (OO1_B, OO1_L, OO1_H));
        var temp_OO2_2 = TranslationWGS84.ConvertGeodeticToLocalTangent(
            (OO2_B, OO2_L, OO2_H),
            (OO1_B, OO1_L, OO1_H));

        /*
         * 
         * ЗАМЕЧАНИЕ
         * Метод GetTargetCoordinates работает при условии что x - направление на Север, y - Направление на Восток
         * 
         */

        (var B2, var T2) = TargetLocator.GetTargetCoordinates(
                            temp_OO1_2.north, temp_OO1_2.east, temp_OO1_2.up,
                            temp_OO2_2.north, temp_OO2_2.east, temp_OO2_2.up,
                            L1, L2, L3, a, b, aa, bb, cc);

        var gB2 = TranslationWGS84.ConvertLocalTangentToGeodetic((B2.Y, B2.X, B2.Z), (OO1_B, OO1_L, OO1_H));
        var gT2 = TranslationWGS84.ConvertLocalTangentToGeodetic((T2.Y, T2.X, T2.Z), (OO1_B, OO1_L, OO1_H));

        Console.WriteLine($"B = {(gB2.B / Math.PI * 180):f7}, {(gB2.L / Math.PI * 180):f7}, {gB2.H:f2}");
        Console.WriteLine($"T = {(gT2.B / Math.PI * 180):f7}, {(gT2.L / Math.PI * 180):f7}, {gT2.H:f2}");
        Console.WriteLine();



        Console.WriteLine();

        return true;
    }
}

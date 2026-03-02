

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


    }

    private static bool test_1((double B, double L, double H) Point)
    {
        var temp = TranslationWGS84.ConvertToGaussKrueger(Point.B, Point.L);
        var point_1 = TranslationWGS84.ConvertFromGaussKrueger(temp.x, temp.y);

        var temp2 = TranslationWGS84.ConvertToGaussKrueger(point_1.B, point_1.L);

        Console.WriteLine(Math.Sqrt(Math.Pow(temp.x - temp2.x,2) + Math.Pow(temp.y - temp2.y, 2)));

        return true;
    }

    private static bool test_2((double B, double L, double H) A, (double B, double L, double H) B)
    {
        var temp_A = TranslationWGS84.ConvertToGaussKrueger(A.B, A.L);
        var temp_B = TranslationWGS84.ConvertToGaussKrueger(B.B, B.L);

        Console.WriteLine(Math.Sqrt(Math.Pow(temp_A.x - temp_B.x, 2) + Math.Pow(temp_A.y - temp_B.y, 2)));

        return true;
    }
}

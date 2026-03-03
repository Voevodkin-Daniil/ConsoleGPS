using System.Text.RegularExpressions;


namespace ConsoleGPS
{
    public class ParserGMS
    {
        /// <summary>
        /// Парсит строку в формате "градусы минуты секунды" в десятичные градусы
        /// Поддерживает форматы: 53°12'05.41", 53°12'5.41", 53 12 5.41
        /// </summary>
        public static double ParseToDecimal(string dmsString)
        {
            // Регулярное выражение для извлечения градусов, минут и секунд
            var pattern = @"(\d+)[°\s](\d+)[′'\s](\d+[.,]?\d*)[\""\s]?";
            var match = Regex.Match(dmsString, pattern);

            if (!match.Success)
                throw new FormatException("Неверный формат координат. Ожидается: 53°12'05.41\"");

            // Парсим значения
            int degrees = int.Parse(match.Groups[1].Value);
            int minutes = int.Parse(match.Groups[2].Value);
            double seconds = double.Parse(match.Groups[3].Value.Replace('.', ','));

            // Конвертируем в десятичные градусы
            double decimalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;

            return Math.Round(decimalDegrees, 6);
        }

        /// <summary>
        /// Конвертирует десятичные градусы в формат "градусы минуты секунды"
        /// </summary>
        public static string ParseToDMS(double decimalDegrees)
        {
            // Определяем знак
            string sign = decimalDegrees >= 0 ? "" : "-";
            decimalDegrees = Math.Abs(decimalDegrees);

            // Вычисляем градусы
            int degrees = (int)Math.Floor(decimalDegrees);

            // Вычисляем минуты
            double minutesDouble = (decimalDegrees - degrees) * 60;
            int minutes = (int)Math.Floor(minutesDouble);

            // Вычисляем секунды
            double seconds = (minutesDouble - minutes) * 60;

            // Форматируем результат
            return $"{sign}{degrees}°{minutes:00}'{seconds:00.00}\"";
        }
    }

    #region УДАЛИТЬ ПОТОМ
    /// <summary>
    /// Предоставляет методы для триангуляции координат Наблюдателя и цели на основе данных с двух наблюдательных пунктов.
    /// </summary>
    internal static class TargetLocator
    {
        /// <summary>
        /// Вычисляет положительный остаток от деления (математический модуль).
        /// </summary>
        /// <param name="x">Делимое</param>
        /// <param name="m">Модуль (делитель)</param>
        /// <returns>Положительное значение остатка от деления</returns>
        private static double Mod(double x, double m)
        {
            double result = x % m;
            return result < 0 ? result + m : result;
        }

        /// <summary>
        /// Вычисляет среднюю точку и расстояние между двумя точками на плоскости.
        /// </summary>
        /// <param name="point1">Координаты первой точки (X, Y, Z)</param>
        /// <param name="point2">Координаты второй точки (X, Y, Z)</param>
        /// <returns>Кортеж, содержащий координаты средней точки и расстояние между исходными точками</returns>
        public static (double X, double Y, double L) AveragePoints((double X, double Y, double Z) point1, (double X, double Y, double Z) point2)
        {
            return
            (
                (point1.X + point2.X) / 2,
                (point1.Y + point2.Y) / 2,
                Math.Sqrt(Math.Pow(point1.X - point2.X, 2) + Math.Pow(point1.Y - point2.Y, 2))
            );
        }

        /// <summary>
        /// Вычисляет среднюю точку и расстояние между двумя точками на плоскости.
        /// </summary>
        /// <param name="point1">Координаты первой точки (X, Y)</param>
        /// <param name="point2">Координаты второй точки (X, Y)</param>
        /// <returns>Кортеж, содержащий координаты средней точки и расстояние между исходными точками</returns>
        public static (double X, double Y, double L) AveragePoints((double X, double Y) point1, (double X, double Y) point2)
            => AveragePoints((point1.X, point1.Y, 0), (point2.X, point2.Y, 0));

        /// <summary>
        /// Находит точки пересечения двух окружностей на плоскости.
        /// </summary>
        /// <param name="x1">X-координата центра первой окружности</param>
        /// <param name="y1">Y-координата центра первой окружности</param>
        /// <param name="x2">X-координата центра второй окружности</param>
        /// <param name="y2">Y-координата центра второй окружности</param>
        /// <param name="radius1">Радиус первой окружности</param>
        /// <param name="radius2">Радиус второй окружности</param>
        /// <param name="height">Z-координата точек пересечения (по умолчанию 0)</param>
        /// <returns>Две точки пересечения окружностей с заданной высотой</returns>
        /// <exception cref="Exception">Возникает, если окружности не пересекаются или совпадают</exception>
        public static ((double X, double Y, double Z) A, (double X, double Y, double Z) B) FindCircleIntersection(
            double x1, double y1, double x2, double y2,
            double radius1, double radius2, double height = 0)
        {
            double dx = x2 - x1;
            double dy = y2 - y1;
            double distance = Math.Sqrt(dx * dx + dy * dy);

            // Проверка геометрической возможности пересечения
            if (distance > radius1 + radius2)
                throw new Exception("Окружности не пересекаются: расстояние между центрами превышает сумму радиусов");

            if (distance < Math.Abs(radius1 - radius2))
                throw new Exception("Окружности не пересекаются: одна окружность находится полностью внутри другой");

            if (distance == 0 && radius1 == radius2)
                throw new Exception("Окружности совпадают: бесконечное множество точек пересечения");

            // Вычисление параметров точек пересечения
            double a = (Math.Pow(radius1, 2) - Math.Pow(radius2, 2) + Math.Pow(distance, 2)) / (2 * distance);
            double h = Math.Sqrt(Math.Pow(radius1, 2) - Math.Pow(a, 2));

            double x_mid = x1 + a * dx / distance;
            double y_mid = y1 + a * dy / distance;

            // Формирование двух точек пересечения
            return (
                (x_mid - h * dy / distance, y_mid + h * dx / distance, height),
                (x_mid + h * dy / distance, y_mid - h * dx / distance, height)
                );
        }

        /// <summary>
        /// Вычисляет координаты Наблюдателя и наземной цели методом триангуляции по данным с двух наблюдателей.
        /// </summary>
        /// <param name="O1_X">Координата X первого опорного оъекта</param>
        /// <param name="O1_Y">Координата Y первого опорного оъекта</param>
        /// <param name="O1_Z">Высота первого опорного оъекта</param>
        /// <param name="O2_X">Координата X второго опорного оъекта</param>
        /// <param name="O2_Y">Координата Y второго опорного оъекта</param>
        /// <param name="O2_Z">Высота второго опорного оъекта</param>
        /// <param name="L1">Наклонная дальность от первого опорного оъекта до Наблюдателя</param>
        /// <param name="L2">Наклонная дальность от второго опорного оъекта до Наблюдателя</param>
        /// <param name="L3">Наклонная дальность от Наблюдателя до цели</param>
        /// <param name="a">Горизонтальный угол между направлением на цель и линией О1-Наблюдатель (в градусах)</param>
        /// <param name="b">Горизонтальный угол между направлением на цель и линией О2-Наблюдатель (в градусах)</param>
        /// <param name="aa">Угол места от первого опорного объекта на Наблюдатель (в градусах)</param>
        /// <param name="bb">Угол места от второго опорного объекта на Наблюдатель (в градусах)</param>
        /// <param name="c">Угол места от Наблюдателя на цель (в градусах)</param>
        /// <returns>
        /// Кортеж из двух точек:
        /// - A: координаты Наблюдателя (X, Y, Z)
        /// - B: координаты цели (X, Y, Z)
        /// </returns>
        public static ((double X, double Y, double Z) A, (double X, double Y, double Z) B) GetTargetCoordinates(
            double O1_X, double O1_Y, double O1_Z, double O2_X, double O2_Y, double O2_Z,
            double L1, double L2, double L3, double a, double b, double aa, double bb, double cc)
        {
            // Расчёт высоты Наблюдателя как среднего значения от двух наблюдателей
            double H = ((O1_Z - Math.Sin(aa * Math.PI / 180) * L1) +
                        (O2_Z - Math.Sin(bb * Math.PI / 180) * L2)) / 2;

            // Горизонтальные проекции дальностей
            double l1 = L1 * Math.Cos(aa * Math.PI / 180);
            double l2 = L2 * Math.Cos(bb * Math.PI / 180);

            // Определение местоположения Наблюдателя (пересечение двух окружностей)
            var BPLA2 = FindCircleIntersection(O1_X, O1_Y, O2_X, O2_Y, l1, l2, 0);
            // Выбор правильного решения на основе разности углов
            var BPLA = (Mod((Math.PI + Math.PI * (a - b) / 180), (2 * Math.PI)) - Math.PI) > 0 ? BPLA2.A : BPLA2.B;

            // Горизонтальная проекция дальности от Наблюдателя до цели
            double l3 = L3 * Math.Cos(cc * Math.PI / 180);

            // Расчёт горизонтальных дальностей от наблюдателей до цели
            double r1 = Math.Sqrt(l1 * l1 + l3 * l3 - 2 * l1 * l3 * Math.Cos(Math.PI * Math.Abs(a) / 180));
            double r2 = Math.Sqrt(l2 * l2 + l3 * l3 - 2 * l2 * l3 * Math.Cos(Math.PI * Math.Abs(b) / 180));

            // Нахождение возможных позиций цели относительно Наблюдателя и опорных объектов
            var points1 = FindCircleIntersection(BPLA.X, BPLA.Y, O1_X, O1_Y, l3, r1);
            var points2 = FindCircleIntersection(BPLA.X, BPLA.Y, O2_X, O2_Y, l3, r2);

            // Поиск наилучшего решения методом усреднения всех комбинаций
            var avgPoints = new List<(double X, double Y, double L)>
                {
                    AveragePoints(points1.A, points2.A),
                    AveragePoints(points1.A, points2.B),
                    AveragePoints(points1.B, points2.A),
                    AveragePoints(points1.B, points2.B)
                };

            // Выбор комбинации с минимальным расстоянием между точками (наиболее согласованное решение)
            avgPoints.Sort((p1, p2) => p1.L.CompareTo(p2.L));
            var bestMatch = avgPoints[0];

            return
            (
                (BPLA.X, BPLA.Y, H), // Координаты Наблюдателя
                (bestMatch.X, bestMatch.Y, H + L3 * Math.Sin(cc * Math.PI / 180)) // Координаты цели
            );
        }

        /// <summary>
        /// Альтернативный метод вычисления координат цели с использованием векторного вращения.
        /// </summary>
        /// <param name="O1_X">Координата X первого опорного оъекта</param>
        /// <param name="O1_Y">Координата Y первого опорного оъекта</param>
        /// <param name="O1_Z">Высота первого опорного оъекта</param>
        /// <param name="O2_X">Координата X второго опорного оъекта</param>
        /// <param name="O2_Y">Координата Y второго опорного оъекта</param>
        /// <param name="O2_Z">Высота второго опорного оъекта</param>
        /// <param name="L1">Наклонная дальность от первого опорного оъекта до Наблюдателя</param>
        /// <param name="L2">Наклонная дальность от второго опорного оъекта до Наблюдателя</param>
        /// <param name="L3">Наклонная дальность от Наблюдателя до цели</param>
        /// <param name="a">Горизонтальный угол между направлением на цель и линией О1-Наблюдатель (в градусах)</param>
        /// <param name="b">Горизонтальный угол между направлением на цель и линией О2-Наблюдатель (в градусах)</param>
        /// <param name="aa">Угол места от первого опорного оъекта на Наблюдатель (в градусах)</param>
        /// <param name="bb">Угол места от второго опорного оъекта на Наблюдатель (в градусах)</param>
        /// <param name="c">Угол места от Наблюдателя на цель (в градусах)</param>
        /// <returns>
        /// Кортеж из двух точек:
        /// - A: координаты Наблюдатель (X, Y, Z)
        /// - B: координаты цели (X, Y, Z)
        /// </returns>
        public static ((double X, double Y, double Z) A, (double X, double Y, double Z) B) GetTargetCoordinatesAlternative(
            double O1_X, double O1_Y, double O1_Z, double O2_X, double O2_Y, double O2_Z,
            double L1, double L2, double L3, double a, double b, double aa, double bb, double c)
        {
            // Расчёт высоты Наблюдателя как среднего значения от двух наблюдателей
            double H = ((O1_Z - Math.Sin(aa * Math.PI / 180) * L1) +
                        (O2_Z - Math.Sin(bb * Math.PI / 180) * L2)) / 2;

            // Горизонтальные проекции дальностей от опорных объектов до Наблюдателя
            double l1 = L1 * Math.Cos(aa * Math.PI / 180);
            double l2 = L2 * Math.Cos(bb * Math.PI / 180);

            // Определение местоположения Наблюдателя (пересечение двух окружностей)
            var BPLA2 = FindCircleIntersection(O1_X, O1_Y, O2_X, O2_Y, l1, l2, 0);
            var BPLA = (Mod((Math.PI + Math.PI * (a - b) / 180), (2 * Math.PI)) - Math.PI) > 0 ? BPLA2.B : BPLA2.A;

            // Горизонтальная проекция дальности от Наблюдателя до цели
            double l3 = L3 * Math.Cos(c * Math.PI / 180);

            // Векторный метод: вращение векторов от Наблюдателя к наблюдателям на заданные углы
            double x1 = O1_X - BPLA.X;
            double y1 = O1_Y - BPLA.Y;
            double x2 = O2_X - BPLA.X;
            double y2 = O2_Y - BPLA.Y;

            // Поворот векторов и масштабирование для получения координат цели
            double rez_X1 = BPLA.X + l3 * (x1 * Math.Cos(Math.PI * -a / 180) - y1 * Math.Sin(Math.PI * -a / 180)) / l1;
            double rez_X2 = BPLA.X + l3 * (x2 * Math.Cos(Math.PI * -b / 180) - y2 * Math.Sin(Math.PI * -b / 180)) / l2;
            double rez_Y1 = BPLA.Y + l3 * (x1 * Math.Sin(Math.PI * -a / 180) + y1 * Math.Cos(Math.PI * -a / 180)) / l1;
            double rez_Y2 = BPLA.Y + l3 * (x2 * Math.Sin(Math.PI * -b / 180) + y2 * Math.Cos(Math.PI * -b / 180)) / l2;

            return
            (
                (BPLA.X, BPLA.Y, H), // Координаты Наблюдатель
                ((rez_X1 + rez_X2) / 2, (rez_Y1 + rez_Y2) / 2, H + L3 * Math.Sin(c * Math.PI / 180)) // Координаты цели (усреднение двух решений)
            );
        }
    }
    #endregion
}

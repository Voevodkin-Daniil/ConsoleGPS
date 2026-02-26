using System;

namespace ConsoleGPS
{
    public enum EncodingDatum
    {
        WGS_84,      // WGS-84
        GSK_2011,    // ГСК-2011
        PZ_90,       // ПЗ-90.11
        SK           // СК-42 / СК-95 (эллипсоид Красовского)
    }

    internal static class TranslationWGS84
    {
        #region Параметры эллипсоидов
        /// <summary>
        /// Возвращает параметры референц-эллипсоида для заданной системы координат
        /// </summary>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Кортеж (a, f), где a - большая полуось в метрах, f - сжатие</returns>
        /// <exception cref="ArgumentException">Выбрасывается при неизвестной системе координат</exception>
        private static (double a, double f) GetEllipsoidParams(EncodingDatum datum)
        {
            return datum switch
            {
                EncodingDatum.WGS_84 => (6378137.0, 1.0 / 298.257223563),
                EncodingDatum.GSK_2011 => (6378136.5, 1.0 / 298.2564151),
                EncodingDatum.PZ_90 => (6378136.0, 1.0 / 298.25784),
                EncodingDatum.SK => (6378245.0, 1.0 / 298.3),
                _ => throw new ArgumentException("Неизвестная система координат")
            };
        }
        #endregion

        #region Преобразование геодезических <-> декартовых координат
        /// <summary>
        /// Преобразует геодезические координаты (широта, долгота, высота) в геоцентрические прямоугольные координаты (X, Y, Z)
        /// </summary>
        /// <param name="latitude">Геодезическая широта в радианах</param>
        /// <param name="longitude">Геодезическая долгота в радианах</param>
        /// <param name="datum">Геодезический датум</param>
        /// <param name="altitude">Высота над эллипсоидом в метрах (по умолчанию 0)</param>
        /// <returns>Массив [X, Y, Z] в метрах</returns>
        public static (double X, double Y, double Z) ConvertGeodeticToCartesian(double latitude, double longitude, EncodingDatum datum, double altitude = 0)
        {
            var (a, f) = GetEllipsoidParams(datum);
            double e2 = 2 * f - f * f; // квадрат первого эксцентриситета
            double sinLat = Math.Sin(latitude);
            double N = a / Math.Sqrt(1 - e2 * sinLat * sinLat);

            double x = (N + altitude) * Math.Cos(latitude) * Math.Cos(longitude);
            double y = (N + altitude) * Math.Cos(latitude) * Math.Sin(longitude);
            double z = ((1 - e2) * N + altitude) * sinLat;

            return (x, y, z );
        }

        /// <summary>
        /// Преобразует геодезические координаты (широта, долгота, высота) в геоцентрические прямоугольные координаты (X, Y, Z)
        /// </summary>
        /// <param name="coordinates">Массив [широта (рад), долгота (рад), высота (м)] (высота опциональна)</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [X, Y, Z] в метрах</returns>
        /// <exception cref="ArgumentException">Выбрасывается, если массив содержит менее 2 элементов</exception>
        public static (double X, double Y, double Z) ConvertGeodeticToCartesian((double B, double L) var, EncodingDatum datum, double altitude = 0)
        {
            return ConvertGeodeticToCartesian(var.B, var.L, datum, altitude);
        }

        /// <summary>
        /// Преобразует геоцентрические прямоугольные координаты (X, Y, Z) в геодезические координаты (широта, долгота, высота)
        /// </summary>
        /// <param name="x">Координата X в метрах</param>
        /// <param name="y">Координата Y в метрах</param>
        /// <param name="z">Координата Z в метрах</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [широта (рад), долгота (рад), высота (м)]</returns>
        public static (double B, double L, double H) ConvertCartesianToGeodetic(double x, double y, double z, EncodingDatum datum)
        {
            var (a, f) = GetEllipsoidParams(datum);
            double e2 = 2 * f - f * f;
            double eps = 1e-10; // допуск итерации

            double p = Math.Sqrt(x * x + y * y); // расстояние до оси Z

            // Случай точки на оси вращения
            if (p < eps)
            {
                double lat2 = (z >= 0) ? Math.PI / 2 : -Math.PI / 2;
                double lon2 = 0;
                double sinLat = Math.Sin(lat2);
                double N = a / Math.Sqrt(1 - e2 * sinLat * sinLat);
                double h2 = z * sinLat - N * (1 - e2 * sinLat * sinLat);
                return (lat2, lon2, h2);
            }

            // Долгота
            double lon = Math.Atan2(y, x);
            if (lon < 0) lon += 2 * Math.PI;

            // Итеративное вычисление широты по ГОСТ (алгоритм 5.1.2)
            double r = Math.Sqrt(p * p + z * z);
            double c = Math.Asin(z / r);
            double p1 = e2 * a / (2 * r);

            double s1 = 0;
            double s2;
            double lat = 0;
            do
            {
                double b = c + s1;
                s2 = Math.Asin(p1 * Math.Sin(2 * b) / Math.Sqrt(1 - e2 * Math.Pow(Math.Sin(b), 2)));
                lat = b;
                if (Math.Abs(s2 - s1) < eps) break;
                s1 = s2;
            } while (true);

            // Высота
            double sinLatFin = Math.Sin(lat);
            double Nfin = a / Math.Sqrt(1 - e2 * sinLatFin * sinLatFin);
            double h = p * Math.Cos(lat) + z * sinLatFin - Nfin * (1 - e2 * sinLatFin * sinLatFin);

            return (lat, lon, h);
        }

        /// <summary>
        /// Преобразует геоцентрические прямоугольные координаты (X, Y, Z) в геодезические координаты (широта, долгота, высота)
        /// </summary>
        /// <param name="cartesian">Массив [X, Y, Z] в метрах</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [широта (рад), долгота (рад), высота (м)]</returns>
        /// <exception cref="ArgumentException">Выбрасывается, если массив содержит менее 3 элементов</exception>
        public static (double B, double L, double H) ConvertCartesianToGeodetic((double X, double Y, double Z) var, EncodingDatum datum)
        {
            return ConvertCartesianToGeodetic(var.X, var.Y, var.Z, datum);
        }
        #endregion

        #region Преобразования между датумами (по ГОСТ 32453-2017)
        /// <summary>
        /// Преобразует геоцентрические координаты из системы ПЗ-90.11 в WGS-84(G1150)
        /// </summary>
        /// <param name="pz90">Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе WGS-84 (метры)</returns>
        public static (double X, double Y, double Z) ConvertPZ90ToWGS84((double X, double Y, double Z) pz90)
        {
            return (
            pz90.X + (+2.041066e-8) * pz90.Y + (+1.716240e-8) * pz90.Z - (-0.003),
            (-2.041066e-8) * pz90.X + pz90.Y + (+1.115071e-8) * pz90.Z - (-0.001),
            (-1.716240e-8) * pz90.X + (-1.115071e-8) * pz90.Y + pz90.Z - (0.000)
            );
        }

        /// <summary>
        /// Преобразует геоцентрические координаты из системы WGS-84(G1150) в ПЗ-90.11
        /// </summary>
        /// <param name="wgs84">Массив [X, Y, Z] в системе WGS-84 (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</returns>
        public static (double X, double Y, double Z) ConvertWGS84ToPZ90((double X, double Y, double Z) wgs84)
        {
            return (
            wgs84.X + (-2.041066e-8) * wgs84.Y + (-1.716240e-8) * wgs84.Z + (-0.003),
            (+2.041066e-8) * wgs84.X + wgs84.Y + (-1.115071e-8) * wgs84.Z + (+0.001),
            (+1.716240e-8) * wgs84.X + (+1.115071e-8) * wgs84.Y + wgs84.Z + (0.000)
            );
        }

        /// <summary>
        /// Преобразует геоцентрические координаты из системы ПЗ-90.11 в ГСК-2011
        /// </summary>
        /// <param name="pz90">Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе ГСК-2011 (метры)</returns>
        public static (double X, double Y, double Z) ConvertPZ90ToGSK2011((double X, double Y, double Z) pz90)
        {
            return (
            pz90.X + (-2.56951e-10) * pz90.Y + (-9.21146e-11) * pz90.Z - (0.000),
            (+2.569513e-10) * pz90.X + pz90.Y + (-2.72465e-9) * pz90.Z - (+0.014),
            (+9.211460e-11) * pz90.X + (-2.72465e-9) * pz90.Y + pz90.Z - (-0.008)
            );
        }

        /// <summary>
        /// Преобразует геоцентрические координаты из системы ГСК-2011 в ПЗ-90.11
        /// </summary>
        /// <param name="gsk2011">Массив [X, Y, Z] в системе ГСК-2011 (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</returns>
        public static (double X, double Y, double Z) ConvertGSK2011ToPZ90((double X, double Y, double Z) gsk2011)
        {
            return (
            gsk2011.X + (+2.569513e-10) * gsk2011.Y + (+9.211460e-11) * gsk2011.Z + (0.000),
            (-2.569513e-10) * gsk2011.X + gsk2011.Y + (+2.724653e-9) * gsk2011.Z + (+0.014),
            (-9.211460e-11) * gsk2011.X + (-2.724653e-9) * gsk2011.Y + gsk2011.Z + (-0.008)
            );
        }

        /// <summary>
        /// Преобразует геоцентрические координаты из системы ПЗ-90.11 в СК-42/95 (эллипсоид Красовского)
        /// </summary>
        /// <param name="pz90">Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе СК (метры)</returns>
        public static (double X, double Y, double Z) ConvertPZ90ToSK((double X, double Y, double Z) pz90)
        {
            return (
            pz90.X + (+6.506684e-7) * pz90.Y + (+1.716240e-8) * pz90.Z - (+24.457),
            (-6.506684e-7) * pz90.X + pz90.Y + (+1.115071e-8) * pz90.Z - (-130.784),
            (-1.716240e-8) * pz90.X + (-1.115071e-8) * pz90.Y + pz90.Z - (-81.538)
            );
        }

        /// <summary>
        /// Преобразует геоцентрические координаты из системы СК-42/95 (эллипсоид Красовского) в ПЗ-90.11
        /// </summary>
        /// <param name="sk">Массив [X, Y, Z] в системе СК (метры)</param>
        /// <returns>Массив [X, Y, Z] в системе ПЗ-90.11 (метры)</returns>
        public static (double X, double Y, double Z) ConvertSKToPZ90((double X, double Y, double Z) sk)
        {
            return (
            sk.X + (-6.506684e-7) * sk.Y + (-1.716240e-8) * sk.Z + (+24.457),
            (+6.506684e-7) * sk.X + sk.Y + (-1.115071e-8) * sk.Z + (-130.784),
            (+1.716240e-8) * sk.X + (+1.115071e-8) * sk.Y + sk.Z + (-81.538)
            );
        }
        #endregion

        #region Проекция Гаусса-Крюгера (эллипсоид Красовского)
        /// <summary>
        /// Преобразует геодезические координаты WGS-84 в плоские прямоугольные координаты проекции Гаусса-Крюгера (зональная система, эллипсоид Красовского)
        /// </summary>
        /// <param name="latitudeWgs84">Геодезическая широта в WGS-84 (радианы)</param>
        /// <param name="longitudeWgs84">Геодезическая долгота в WGS-84 (радианы)</param>
        /// <returns>Массив [x, y] в метрах, где x - северное смещение (ордината), y - восточное смещение с номером зоны</returns>
        public static (double x, double y) ConvertToGaussKrueger(double latitudeWgs84, double longitudeWgs84)
        {
            // Цепочка: WGS84 (геод.) -> декартовы WGS84 -> PZ-90 -> SK -> геод. SK
            var wgsCart = ConvertGeodeticToCartesian(latitudeWgs84, longitudeWgs84, EncodingDatum.WGS_84);
            var pzCart = ConvertWGS84ToPZ90(wgsCart);
            var skCart = ConvertPZ90ToSK(pzCart);
            var skGeo = ConvertCartesianToGeodetic(skCart, EncodingDatum.SK);
            double B = skGeo.B; // широта в радианах (эллипсоид Красовского)
            double L = skGeo.L; // долгота в радианах

            // Прямая проекция Гаусса-Крюгера (формулы 24-26 ГОСТ)
            double L_deg = L * 180.0 / Math.PI;
            int n = (int)((6 + L_deg) / 6);
            double L0 = (3 + 6 * (n - 1)) * Math.PI / 180.0; // осевой меридиан в радианах
            double l = L - L0;

            double sinB = Math.Sin(B);
            double sin2B = Math.Sin(2 * B);
            double sinB2 = sinB * sinB;
            double sinB4 = sinB2 * sinB2;
            double sinB6 = sinB4 * sinB2;
            double l2 = l * l;

            // Вычисление x по формуле (24)
            double term1 = 16002.8900 + 66.9607 * sinB2 + 0.3515 * sinB4;
            double term2 = 1594561.25 + 5336.535 * sinB2 + 26.790 * sinB4 + 0.149 * sinB6;
            double term3 = 672483.4 - 811219.9 * sinB2 + 5420.0 * sinB4 - 10.6 * sinB6;
            double term4 = 278194 - 830174 * sinB2 + 572434 * sinB4 - 16010 * sinB6;
            double term5 = 109500 - 574700 * sinB2 + 863700 * sinB4 - 398600 * sinB6;

            double x = 6367558.4968 * B
                - sin2B * (term1
                    - l2 * (term2
                        + l2 * (term3
                            + l2 * (term4
                                + l2 * term5))));

            // Вычисление y по формуле (25)
            double y = (5 + 10 * n) * 100000
                + l * Math.Cos(B) * (6378245 + 21346.1415 * sinB2 + 107.1590 * sinB4 + 0.5977 * sinB6
                + l2 * (1070204.16 - 2136826.66 * sinB2 + 17.98 * sinB4 - 11.99 * sinB6
                + l2 * (270806 - 1523417 * sinB2 + 1327645 * sinB4 - 21701 * sinB6
                + l2 * (79690 - 866190 * sinB2 + 1730360 * sinB4 - 945460 * sinB6))));

            return ( x, y );
        }

        /// <summary>
        /// Преобразует плоские прямоугольные координаты проекции Гаусса-Крюгера (зональная система, эллипсоид Красовского) в геодезические координаты WGS-84
        /// </summary>
        /// <param name="x">Северное смещение (ордината) в метрах</param>
        /// <param name="y">Восточное смещение с номером зоны в метрах (первые цифры - номер зоны)</param>
        /// <returns>Массив [широта (рад), долгота (рад), высота (м)] в системе WGS-84</returns>
        public static (double B, double L, double H) ConvertFromGaussKrueger(double x, double y)
        {
            // Обратная проекция (формулы 29-36 ГОСТ)
            int n = (int)(y * 1e-6);
            double beta = x / 6367558.4968;

            // Вычисление B0 по формуле (32)
            double sinBeta = Math.Sin(beta);
            double sinBeta2 = sinBeta * sinBeta;
            double sinBeta4 = sinBeta2 * sinBeta2;
            double sin2Beta = Math.Sin(2 * beta);

            double B0 = beta + sin2Beta * (0.00252588685 - 0.00001491860 * sinBeta2 + 0.00000011904 * sinBeta4);

            // Теперь используем B0 для вычисления ΔB и l
            double sinB0 = Math.Sin(B0);
            double sinB02 = sinB0 * sinB0;
            double sinB04 = sinB02 * sinB02;
            double sinB06 = sinB04 * sinB02;
            double sin2B0 = Math.Sin(2 * B0);
            double cosB0 = Math.Cos(B0);

            double z0 = (y - (10 * n + 5) * 100000) / (6378245 * cosB0);
            double z02 = z0 * z0;

            // Формула (33) для ΔB
            double dB = -z02 * sin2B0 * (0.251684631 - 0.003369263 * sinB02 + 0.00001127 * sinB04
                - z02 * (0.10500614 - 0.04559916 * sinB02 + 0.00228901 * sinB04 - 0.00002987 * sinB06
                - z02 * (0.042858 - 0.025318 * sinB02 + 0.014346 * sinB04 - 0.001264 * sinB06
                - z02 * (0.01672 - 0.00630 * sinB02 + 0.01188 * sinB04 - 0.00328 * sinB06))));

            // Формула (34) для l
            double l = z0 * (1 - 0.0033467108 * sinB02 - 0.0000056002 * sinB04 - 0.0000000187 * sinB06
                - z02 * (0.16778975 + 0.16273586 * sinB02 - 0.00052490 * sinB04 - 0.00000846 * sinB06
                - z02 * (0.0420025 + 0.1487407 * sinB02 + 0.0059420 * sinB04 - 0.0000150 * sinB06
                - z02 * (0.01225 + 0.09477 * sinB02 + 0.03282 * sinB04 - 0.00034 * sinB06
                - z02 * (0.0038 + 0.0524 * sinB02 + 0.0482 * sinB04 - 0.0032 * sinB06)))));

            // Геодезические координаты в системе СК (эллипсоид Красовского)
            double B_sk = B0 + dB;
            double L_sk = (6 * (n - 0.5)) * Math.PI / 180.0 + l;

            // Обратная цепочка: СК (геод.) -> декартовы СК -> PZ-90 -> WGS84 -> геод. WGS84
            var skCart = ConvertGeodeticToCartesian(B_sk, L_sk, EncodingDatum.SK);
            var pzCart = ConvertSKToPZ90(skCart);
            var wgsCart = ConvertPZ90ToWGS84(pzCart);
            var wgsGeo = ConvertCartesianToGeodetic(wgsCart, EncodingDatum.WGS_84);

            return wgsGeo; // [широта (рад), долгота (рад), высота (м)]
        }
        #endregion
    }
}
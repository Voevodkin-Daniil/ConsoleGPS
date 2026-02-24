using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Intrinsics.Arm;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleGPS
{

    /*
     * ПЕРЕПИСАТЬ ВЛГОРИТМ И СДЕЛАТЬ ЕГО КОРЕКТНЫМ
     */


    internal static class TestLib
    {

        /// <summary>
        /// Преобразует географические координаты (широта, долгота, высота) в геоцентрические координаты XYZ
        /// </summary>
        /// <param name="latitude">Широта в радианах</param>
        /// <param name="longitude">Долгота в радианах</param>
        /// <param name="datum">Геодезический датум (система координат)</param>
        /// <param name="altitude">Высота над эллипсоидом в метрах (по умолчанию 0)</param>
        /// <returns>Массив [X, Y, Z] геоцентрических координат в метрах</returns>
        public static double[] ConvertGeodeticToCartesian(double latitude, double longitude, EncodingDatum datum, double altitude = 0)
        {
            // Константы эллипсоида для разных датумов
            double semiMajorAxis = datum switch
            {
                EncodingDatum.WGS_84 => 6378137.0,
                EncodingDatum.GSK_2011 => 6378136.5,
                EncodingDatum.PZ_90 => 6378136,
                EncodingDatum.SK => 6378245,
                _ => -1
            };

            double flattening = datum switch
            {
                EncodingDatum.WGS_84 => 1.0 / 298.257223563,
                EncodingDatum.GSK_2011 => 1.0 / 298.257223563, //!!!!!!!!!!!!!!!!!!!!
                EncodingDatum.PZ_90 => 1.0 / 298.25784,
                EncodingDatum.SK => 1.0 / 298.3,
                _ => -1
            };

            double x = 0, y = 0, z = 0;

            // Первый эксцентриситет
            double eccentricity = Math.Sqrt(2 * flattening - flattening * flattening);

            // Радиус кривизны первого вертикала
            double primeVerticalRadius = semiMajorAxis / Math.Sqrt(1 - Math.Pow(eccentricity * Math.Sin(latitude), 2));

            // Преобразование в геоцентрические координаты
            x = (primeVerticalRadius + altitude) * Math.Cos(latitude) * Math.Cos(longitude);
            y = (primeVerticalRadius + altitude) * Math.Cos(latitude) * Math.Sin(longitude);
            z = ((1 - eccentricity * eccentricity) * primeVerticalRadius + altitude) * Math.Sin(latitude);

            return new double[] { x, y, z };
        }

        /// <summary>
        /// Преобразует массив географических координат в геоцентрические координаты XYZ
        /// </summary>
        /// <param name="coordinates">Массив координат [широта, долгота] или [широта, долгота, высота]</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [X, Y, Z] геоцентрических координат</returns>
        public static double[] ConvertGeodeticToCartesian(double[] coordinates, EncodingDatum datum)
        {
            if (coordinates.Length == 2)
                return ConvertGeodeticToCartesian(coordinates[0], coordinates[1], datum);
            return ConvertGeodeticToCartesian(coordinates[0], coordinates[1], datum, coordinates[2]);
        }

        /// <summary>
        /// Преобразует геоцентрические координаты XYZ в географические (широта, долгота, высота)
        /// </summary>
        /// <param name="x">Координата X в метрах</param>
        /// <param name="y">Координата Y в метрах</param>
        /// <param name="z">Координата Z в метрах</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [широта (рад), долгота (рад), высота (м)]</returns>
        public static double[] ConvertCartesianToGeodetic(double x, double y, double z, EncodingDatum datum)
        {
            // Константы эллипсоида
            double semiMajorAxis = datum switch
            {
                EncodingDatum.WGS_84 => 6378137.0,
                EncodingDatum.GSK_2011 => 6378136.5,
                EncodingDatum.PZ_90 => 6378136,
                _ => -1
            };

            double flattening = datum switch
            {
                EncodingDatum.WGS_84 => 1.0 / 298.257223563,
                EncodingDatum.GSK_2011 => 0,
                EncodingDatum.PZ_90 => 1.0 / 298.25784,
                _ => -1
            };

            double latitude = 0, longitude = 0, altitude = 0;

            // Первый эксцентриситет
            double eccentricity = Math.Sqrt(2 * flattening - flattening * flattening);

            // Преобразование
            double horizontalDistance = Math.Sqrt(x * x + y * y);

            if (horizontalDistance == 0) // Точка на оси Z
            {
                latitude = Math.Sign(z) * Math.PI / 2;
                longitude = 0;
                altitude = z * Math.Sin(latitude) - semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity * Math.Sin(latitude), 2));
            }
            else
            {
                // Вычисление долготы
                double longitudeAngle = Math.Abs(Math.Asin(y / horizontalDistance));

                if (y < 0 && x > 0) longitude = 2 * Math.PI - longitudeAngle;
                else if (y < 0 && x < 0) longitude = Math.PI + longitudeAngle;
                else if (y > 0 && x < 0) longitude = Math.PI - longitudeAngle;
                else if (y > 0 && x > 0) longitude = longitudeAngle;
                else if (y == 0 && x > 0) longitude = 0;
                else longitude = Math.PI;

                if (z == 0) // Точка в экваториальной плоскости
                {
                    latitude = 0;
                    altitude = horizontalDistance - semiMajorAxis;
                }
                else
                {
                    // Итеративное вычисление широты
                    double distance = Math.Sqrt(x * x + y * y + z * z);
                    double centralAngle = Math.Asin(z / distance);
                    double parameter = eccentricity * eccentricity * semiMajorAxis / (2 * distance);

                    double previousDelta = 0, currentDelta = 0;
                    double currentLatitude = 0;

                    double tolerance = 0.000001;
                    double delta = 1;

                    do
                    {
                        currentLatitude = centralAngle + previousDelta;
                        currentDelta = Math.Asin(parameter * Math.Sin(2 * currentLatitude) /
                                                  Math.Sqrt(1 - Math.Pow(eccentricity * Math.Sin(currentLatitude), 2)));
                        delta = Math.Abs(currentDelta - previousDelta);
                        previousDelta = currentDelta;
                    } while (delta > tolerance);

                    latitude = currentLatitude;
                    altitude = horizontalDistance * Math.Cos(latitude) + z * Math.Sin(latitude) -
                               semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity * Math.Sin(latitude), 2));
                }
            }

            return new double[] { latitude, longitude, altitude };
        }

        /// <summary>
        /// Преобразует массив геоцентрических координат в географические
        /// </summary>
        /// <param name="cartesianCoordinates">Массив [X, Y, Z] в метрах</param>
        /// <param name="datum">Геодезический датум</param>
        /// <returns>Массив [широта (рад), долгота (рад), высота (м)]</returns>
        public static double[] ConvertCartesianToGeodetic(double[] cartesianCoordinates, EncodingDatum datum)
        {
            return ConvertCartesianToGeodetic(cartesianCoordinates[0], cartesianCoordinates[1], cartesianCoordinates[2], datum);
        }

        /// <summary>
        /// Преобразует координаты из датума WGS84 в датум PZ-90 (геоцентрические координаты)
        /// </summary>
        /// <param name="wgs84Coordinates">Массив [X, Y, Z] в WGS84</param>
        /// <returns>Массив [X, Y, Z] в PZ-90</returns>
        public static double[] ConvertWGS84ToPZ90(double[] wgs84Coordinates)
        {
            return ConvertWGS84ToPZ90(wgs84Coordinates[0], wgs84Coordinates[1], wgs84Coordinates[2]);
        }

        /// <summary>
        /// Преобразует координаты из датума WGS84 в датум PZ-90 (геоцентрические координаты)
        /// </summary>
        /// <param name="xWgs84">X координата в WGS84 (м)</param>
        /// <param name="yWgs84">Y координата в WGS84 (м)</param>
        /// <param name="zWgs84">Z координата в WGS84 (м)</param>
        /// <returns>Массив [X, Y, Z] в PZ-90</returns>
        public static double[] ConvertWGS84ToPZ90(double xWgs84, double yWgs84, double zWgs84)
        {
            double scaleFactor = 1 + (-0.008) * Math.Pow(10, -6);

            double xPz90 = scaleFactor * (
                xWgs84 +
                -2.041066e-8 * yWgs84 +
                -1.716240e-8 * zWgs84
            ) + (-0.003);

            double yPz90 = scaleFactor * (
                -2.041066e-8 * xWgs84 +
                yWgs84 +
                -1.115071e-8 * zWgs84
            ) + (-0.001);

            double zPz90 = scaleFactor * (
                -1.716240e-8 * xWgs84 +
                1.115071e-8 * yWgs84 +
                zWgs84
            ) + 0.0;

            return new double[] { xPz90, yPz90, zPz90 };
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум GSK-2011 (геоцентрические координаты)
        /// </summary>
        /// <param name="pz90Coordinates">Массив [X, Y, Z] в PZ-90</param>
        /// <returns>Массив [X, Y, Z] в GSK-2011</returns>
        public static double[] ConvertPZ90ToGSK2011(double[] pz90Coordinates)
        {
            return ConvertPZ90ToGSK2011(pz90Coordinates[0], pz90Coordinates[1], pz90Coordinates[2]);
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум GSK-2011 (геоцентрические координаты)
        /// </summary>
        /// <param name="xPz90">X координата в PZ-90 (м)</param>
        /// <param name="yPz90">Y координата в PZ-90 (м)</param>
        /// <param name="zPz90">Z координата в PZ-90 (м)</param>
        /// <returns>Массив [X, Y, Z] в GSK-2011</returns>
        public static double[] ConvertPZ90ToGSK2011(double xPz90, double yPz90, double zPz90)
        {
            double scaleFactor = 1 - (-0.006) * Math.Pow(10, -6);

            double xGsk2011 = scaleFactor * (
                xPz90 +
                -2.56951e-10 * yPz90 +
                -9.21146e-11 * zPz90
            ) + 0.000; // Обратите внимание: было -(-0.000), что равно +0.000

            double yGsk2011 = scaleFactor * (
                2.569513e-10 * xPz90 +
                yPz90 +
                -2.72465e-9 * zPz90
            ) - 0.014;

            double zGsk2011 = scaleFactor * (
                9.211460e-11 * xPz90 +
                -2.72465e-9 * yPz90 +
                zPz90
            ) + 0.008; // Обратите внимание: было -(-0.008), что равно +0.008

            return new double[] { xGsk2011, yGsk2011, zGsk2011 };
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум СК (геоцентрические координаты)
        /// </summary>
        /// <param name="pz90Coordinates">Массив [X, Y, Z] в PZ-90</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertPZ90ToSK(double[] pz90Coordinates)
        {
            return ConvertPZ90ToSK(pz90Coordinates[0], pz90Coordinates[1], pz90Coordinates[2]);
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум СК (геоцентрические координаты)
        /// </summary>
        /// <param name="xPz90">X координата в PZ-90 (м)</param>
        /// <param name="yPz90">Y координата в PZ-90 (м)</param>
        /// <param name="zPz90">Z координата в PZ-90 (м)</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertPZ90ToSK(double xPz90, double yPz90, double zPz90)
        {
            // Примечание: не уточнено, СК-45 или СК-92
            double scaleFactor = 1 - (-0.228) * Math.Pow(10, -6);

            double xSk = scaleFactor * (
                xPz90 +
                6.506684e-7 * yPz90 +
                1.716240e-8 * zPz90
            ) - 24.457;

            double ySk = scaleFactor * (
                -6.506684e-7 * xPz90 +
                yPz90 +
                1.115071e-8 * zPz90
            ) + 130.784; // Обратите внимание: было -(-130.784), что равно +130.784

            double zSk = scaleFactor * (
                -1.716240e-8 * xPz90 +
                -1.115071e-8 * yPz90 +
                zPz90
            ) + 81.538; // Обратите внимание: было -(-81.538), что равно +81.538

            return new double[] { xSk, ySk, zSk };
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум WGS84 (геоцентрические координаты)
        /// </summary>
        /// <param name="pz90Coordinates">Массив [X, Y, Z] в PZ-90</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertPZ90ToWGS84(double[] wgs84Coordinates)
        {
            return ConvertPZ90ToWGS84(wgs84Coordinates[0], wgs84Coordinates[1], wgs84Coordinates[2]);
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум WGS84 (геоцентрические координаты)
        /// </summary>
        /// <param name="xPz90">X координата в PZ-90 (м)</param>
        /// <param name="yPz90">Y координата в PZ-90 (м)</param>
        /// <param name="zPz90">Z координата в PZ-90 (м)</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertPZ90ToWGS84(double xPz90, double yPz90, double zPz90)
        {
            // Примечание: не уточнено, СК-45 или СК-92
            double scaleFactor = 1 - (-0.008) * Math.Pow(10, -6);

            double xWgs84 = scaleFactor * (
                xPz90 +
                2.041066e-8 * yPz90 +
                1.716240e-8 * zPz90
            ) + 0.003;

            double yWgs84 = scaleFactor * (
                -2.041066e-8 * xPz90 +
                yPz90 +
                1.115071e-8 * zPz90
            ) + 0.001; // Обратите внимание: было -(-130.784), что равно +130.784

            double zWgs84 = scaleFactor * (
                -1.716240e-8 * xPz90 +
                -1.115071e-8 * yPz90 +
                zPz90
            ) + 0;

            return new double[] { xWgs84, yWgs84, zWgs84 };
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум СК (геоцентрические координаты)
        /// </summary>
        /// <param name="pz90Coordinates">Массив [X, Y, Z] в PZ-90</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertSKToPZ90(double[] skCoordinates)
        {
            return ConvertSKToPZ90(skCoordinates[0], skCoordinates[1], skCoordinates[2]);
        }

        /// <summary>
        /// Преобразует координаты из датума PZ-90 в датум СК (геоцентрические координаты)
        /// </summary>
        /// <param name="xPz90">X координата в PZ-90 (м)</param>
        /// <param name="yPz90">Y координата в PZ-90 (м)</param>
        /// <param name="zPz90">Z координата в PZ-90 (м)</param>
        /// <returns>Массив [X, Y, Z] в СК</returns>
        public static double[] ConvertSKToPZ90(double xSK, double ySK, double zSK)
        {
            // Примечание: не уточнено, СК-45 или СК-95
            double scaleFactor = 1 + (-0.228) * Math.Pow(10, -6);

            double xSk = scaleFactor * (
                xSK +
                -6.506684e-7 * ySK +
                -1.716240e-8 * zSK
            ) + 24.457;

            double ySk = scaleFactor * (
                6.506684e-7 * xSK +
                ySK +
                -1.115071e-8 * zSK
            ) - 130.784; // Обратите внимание: было -(-130.784), что равно +130.784

            double zSk = scaleFactor * (
                1.716240e-8 * xSK +
                1.115071e-8 * ySK +
                zSK
            ) - 81.538; // Обратите внимание: было -(-81.538), что равно +81.538

            return new double[] { xSk, ySk, zSk };
        }

        /// <summary>
        /// Преобразует географические координаты из WGS84 в проекцию Гаусса-Крюгера (плоские координаты)
        /// через промежуточное преобразование в датум Красовского (СК)
        /// </summary>
        /// <param name="latitudeWgs84">Широта в WGS84 (радианы)</param>
        /// <param name="longitudeWgs84">Долгота в WGS84 (радианы)</param>
        /// <returns>Массив [X, Y] плоских координат в проекции Гаусса-Крюгера</returns>
        public static double[] ConvertToGaussKrueger(double latitudeWgs84, double longitudeWgs84)
        {
            // Преобразование из WGS84 в датум Красовского через PZ-90 и СК
            double[] wgs84Cartesian = ConvertGeodeticToCartesian(latitudeWgs84, longitudeWgs84, EncodingDatum.WGS_84);
            double[] pz90Cartesian = ConvertWGS84ToPZ90(wgs84Cartesian);
            double[] skCartesian = ConvertPZ90ToSK(pz90Cartesian);
            double[] skGeodetic = ConvertCartesianToGeodetic(skCartesian, EncodingDatum.SK);

            double latitudeSk = skGeodetic[0];
            double longitudeSk = skGeodetic[1];

            // Проекция Гаусса-Крюгера для эллипсоида Красовского
            double[] result = new double[2];

            // Переводим долготу в градусы для определения зоны
            double longitudeDegrees = longitudeSk * 180 / Math.PI;

            // Номер 6-градусной зоны
            int zoneNumber = (int)((6 + longitudeDegrees) / 6);

            // Расстояние от осевого меридиана зоны (в радианах)
            double centralMeridian = (3 + 6 * (zoneNumber - 1)) / 57.29577951; // перевод в радианы
            double deltaLongitude = longitudeSk - centralMeridian;

            // Предварительные вычисления для оптимизации
            double sinLatitude = Math.Sin(latitudeSk);
            double sinLatitude2 = sinLatitude * sinLatitude;
            double sinLatitude4 = sinLatitude2 * sinLatitude2;
            double sinLatitude6 = sinLatitude4 * sinLatitude2;
            double sin2Latitude = Math.Sin(2 * latitudeSk);
            double deltaLon2 = deltaLongitude * deltaLongitude;

            // Вычисление X (северное смещение) по формуле из ГОСТ Р 32453-2017
            result[0] = 6367558.4968 * latitudeSk -
                        sin2Latitude * (16002.8900 + 66.9607 * sinLatitude2 + 0.3515 * sinLatitude4 -
                        deltaLon2 * (1594561.25 + 5336.535 * sinLatitude2 + 26.790 * sinLatitude4 + 0.149 * sinLatitude6 +
                        deltaLon2 * (672483.4 - 811219.9 * sinLatitude4 + 5420.0 * sinLatitude4 - 10.6 * sinLatitude6 +
                        deltaLon2 * (278194 - 830174 * sinLatitude2 + 572434 * sinLatitude4 - 16010 * sinLatitude6 +
                        deltaLon2 * (109500 - 574700 * sinLatitude2 + 863700 * sinLatitude4 - 398600 * sinLatitude6)))));

            // Вычисление Y (восточное смещение)
            result[1] = (5 + 10 * zoneNumber) * 100000 +
                        deltaLongitude * Math.Cos(latitudeSk) *
                        (6378245 + 21346.1415 * sinLatitude2 + 107.1590 * sinLatitude4 + 0.5977 * sinLatitude6 +
                        deltaLon2 * (1070204.16 - 2136826.66 * sinLatitude2 + 17.98 * sinLatitude4 - 11.99 * sinLatitude6 +
                        deltaLon2 * (270806 - 1523417 * sinLatitude2 + 1327645 * sinLatitude4 - 21701 * sinLatitude6 +
                        deltaLon2 * (79690 - 866190 * sinLatitude2 + 1730360 * sinLatitude4 - 945460 * sinLatitude6))));

            return result;
        }

        /// <summary>
        /// Преобразует проекцию Гаусса-Крюгера в географические координаты из WGS84
        /// </summary>
        /// <param name="latitudeWgs84"></param>
        /// <param name="longitudeWgs84"></param>
        /// <returns></returns>
        public static double[] ConvertFromGaussKrueger(double xSK, double ySK)
        {
            // Преобразование из WGS84 в датум Красовского через PZ-90 и СК


            // Географические координаты стандарт Красавского
            double[] result = new double[2];

            double n = (int)(ySK * Math.Pow(10, -6));
            double beta = xSK / 6367558.4968;

            double sinBeta = Math.Sin(beta);
            double sinBeta2 = sinBeta * sinBeta;
            double sinBeta4 = sinBeta2 * sinBeta2;
            double sinBeta6 = sinBeta4 * sinBeta2;
            double sin2Beta = Math.Sin(2 * beta);

            double B0 = beta + sin2Beta * (0.00252588685 - 0.00001491860 * sinBeta2 + 0.00000011904 * sinBeta4);
            double z0 = (ySK - (10 * n + 5) * 100000) / (6378245 * Math.Cos(B0));
            double z0_2 = z0 * z0;



            double dB = -z0_2 * sin2Beta * (0.251684631 - 0.003369263 * sinBeta2 + 0.00001127 * sinBeta4 -
                         z0_2 * (0.10500614 - 0.04559916 * sinBeta2 + 0.00228901 * sinBeta4 - 0.00002987 * sinBeta6 -
                         z0_2 * (0.042858 - 0.025318 * sinBeta2 + 0.014346 * sinBeta4 - 0.001264 * sinBeta6 -
                         z0_2 * (0.01672 - 0.00630 * sinBeta2 + 0.01188 * sinBeta4 - -0.00328 * sinBeta6))));


            double l = z0 * (1 - 0.0033467108 * sinBeta2 - 0.0000056002 * sinBeta4 - 0.0000000187 * sinBeta6 -
                       z0_2 * (0.16778975 + 0.16273586 * sinBeta2 - 0.00052490 * sinBeta4 - 0.00000846 * sinBeta6 -
                       z0_2 * (0.0420025 + 0.1487407 * sinBeta2 + 0.0059420 * sinBeta4 - 0.0000150 * sinBeta6 -
                       z0_2 * (0.01225 + 0.09477 * sinBeta2 + 0.03282 * sinBeta4 - 0.00034 * sinBeta6 -
                       z0_2 * (0.0038 + 0.0524 * sinBeta2 + 0.0482 * sinBeta4 - 0.0032 * sinBeta6)))));




            result[0] = B0 + dB;
            result[1] = 6 * (n - 0.5) / 57.29577951 + l;

            // Преобразованиие из системы SK в WGS84

            double[] skCartesian = ConvertGeodeticToCartesian(result, EncodingDatum.SK);
            double[] pz90Cartesian = ConvertSKToPZ90(skCartesian);
            double[] wgs84Cartesian = ConvertPZ90ToWGS84(pz90Cartesian);
            double[] wgsGeodetic = ConvertCartesianToGeodetic(wgs84Cartesian, EncodingDatum.WGS_84);

            return wgsGeodetic;
        }
    }
}
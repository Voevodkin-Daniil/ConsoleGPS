using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace ConsoleGPS
{
    internal class ParserGMS
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
}

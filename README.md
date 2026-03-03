# ConsoleGPS — Библиотека и консольная утилита для геодезических вычислений

Проект представляет собой консольное приложение и библиотеку на C# для широкого спектра геодезических расчетов: преобразования координат между различными системами (датумами), проекциями, форматами. В основе лежат алгоритмы, соответствующие стандартам (в частности, ГОСТ 32453-2017).

## Ссылка на репозиторий с оконным приложением для тестирования этой библиотеки
[*Реализация в оконном приложении (WPF)*](https://github.com/Voevodkin-Daniil/WindowedCoordinateDetectionWGS84/ "Связанный проект")

## Источники
- **ГОСТ 32453-2017** — [racurs.ru](https://racurs.ru/downloads/documentation/gost_r_32453-2017.pdf?ysclid=m74cniwn8633815025)
- **Google Earth** — для проверки координат
- **Геодезический калькулятор** — [geoproj.ru](https://geoproj.ru/)
- **Local Tangent Plane Coordinates** — [PSAS Archive](https://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf)

## Основные возможности

### 1. Парсинг и форматирование координат (`ParserGMS`)
*   **`ParseToDecimal()`** — преобразование из формата «Градусы Минуты Секунды» (например, `53°12'05.41"`, `53 12 5.41`) в десятичные градусы.
*   **`ParseToDMS()`** — обратное преобразование из десятичных градусов в строку формата DMS.

### 2. Пересчет между геодезическими и декартовыми координатами (`TranslationWGS84`)
*   Прямое и обратное преобразование для эллипсоидов:
    *   **WGS-84** (Всемирная геодезическая система)
    *   **ГСК-2011** (Государственная система координат РФ)
    *   **ПЗ-90.11** (Параметры Земли 1990 г.)
    *   **СК-42/95** (Система координат 1942/1995 г. на эллипсоиде Красовского)

### 3. Пересчет между различными геодезическими датумами (ГОСТ 32453-2017)
*   WGS-84 (G1150) <=> ПЗ-90.11
*   ПЗ-90.11 <=> ГСК-2011
*   ПЗ-90.11 <=> СК-42/95

### 4. Проекция Гаусса-Крюгера
*   **Прямое преобразование:** WGS-84 (геод.) → Прямоугольные координаты (зональная система СК-42/95, эллипсоид Красовского).
*   **Обратное преобразование:** Прямоугольные координаты (СК-42/95) → Геодезические координаты WGS-84.

### 5. Локальная касательная плоскость (ENU)
*   **`ConvertGeodeticToLocalTangent`** — перевод геодезических координат (WGS-84) в локальную топоцентрическую систему координат «Восток-Север-Вверх» (East-North-Up) относительно заданной опорной точки.
*   **`ConvertLocalTangentToGeodetic`** — обратное преобразование.

## Алгоритмическая схема преобразований датумов
![Image Alg](image/Algoritm.png)

## Структура проекта
*   **`Program.cs`** — Демонстрационный пример. Показывает:
    *   Преобразование координат в проекцию Гаусса-Крюгера и обратно.
    *   Расчет расстояний в локальной касательной плоскости.
    *   Триангуляцию цели на основе реальных данных.
*   **`ParserGMS.cs`** — Парсер форматов DMS и класс `TargetLocator` для триангуляции.
*   **`TranslationWGS84.cs`** — Основной класс геодезических преобразований.

## Пример использования (триангуляция)

```csharp
Console.WriteLine("Преобразование координат демонстрационный вариант");

(double B, double L, double H) OO1 = (53.2015043 * Math.PI / 180.0, 50.1133423 * Math.PI / 180.0, 0);

(double B, double L, double H) OO2 = (53.2000976 * Math.PI / 180.0, 50.1125830 * Math.PI / 180.0, 0);

// Сравнение насколько будет разница при двойном определении
test_1(OO1);

// Проерка что выдаст растояние
test_2(OO1, OO2);
test_3(OO2, OO1);

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

```

## Реализация преобразований геодезических координат (TranslationWGS84)

### Общее описание
Класс `TranslationWGS84` предоставляет комплексный набор методов для преобразования координат между различными геодезическими системами, включая WGS-84, ГСК-2011, ПЗ-90.11 и СК-42/95 (эллипсоид Красовского). Реализация основана на ГОСТ 32453-2017.

### Поддерживаемые системы координат (датумы)

| Датум     | Описание                               | Параметры эллипсоида                 |
|-----------|----------------------------------------|--------------------------------------|
| `WGS_84`  | Всемирная геодезическая система 1984   | a = 6378137.0 м, f = 1/298.257223563 |
| `GSK_2011`| Государственная система координат 2011 | a = 6378136.5 м, f = 1/298.2564151   |
| `PZ_90`   | Параметры Земли 1990.11                | a = 6378136.0 м, f = 1/298.25784     |
| `SK`      | СК-42/95 (эллипсоид Красовского)       | a = 6378245.0 м, f = 1/298.3         |

### 1. Преобразование геодезических и декартовых координат

#### 1.1 Из геодезических в декартовы (Геоцентрические прямоугольные координаты)

```csharp
// Преобразование по отдельным параметрам
public static (double X, double Y, double Z) ConvertGeodeticToCartesian(
    double latitude,      // широта в радианах
    double longitude,     // долгота в радианах  
    EncodingDatum datum,  // датум
    double altitude = 0)  // высота над эллипсоидом в метрах

// Преобразование с использованием кортежа
public static (double X, double Y, double Z) ConvertGeodeticToCartesian(
    (double B, double L) geodetic,  // широта и долгота в радианах
    EncodingDatum datum,             // датум
    double altitude = 0)             // высота над эллипсоидом в метрах
```

**Алгоритм:**
1. Получение параметров эллипсоида для заданного датума
2. Вычисление квадрата первого эксцентриситета: e² = 2f - f²
3. Вычисление радиуса кривизны первого вертикала: N = a / ✓(1 - e²·sin²φ)
4. Преобразование по формулам:
   - X = (N + h)·cos φ·cos λ
   - Y = (N + h)·cos φ·sin λ
   - Z = ((1 - e²)·N + h)·sin φ

#### 1.2 Из декартовых в геодезические

```csharp
// Преобразование по отдельным координатам
public static (double B, double L, double H) ConvertCartesianToGeodetic(
    double x, double y, double z,  // координаты в метрах
    EncodingDatum datum)            // датум

// Преобразование с использованием кортежа
public static (double B, double L, double H) ConvertCartesianToGeodetic(
    (double X, double Y, double Z) cartesian,  // координаты в метрах
    EncodingDatum datum)                        // датум
```

**Алгоритм (итеративный метод по ГОСТ 32453-2017, п. 5.1.2):**
1. Обработка частного случая точки на оси вращения (p < ε)
2. Вычисление долготы: λ = atan2(y, x)
3. Итеративное вычисление широты с использованием вспомогательных величин:
   - p = ✓(x² + y²)
   - r = ✓(p² + z²)
   - c = asin(z / r)
   - p₁ = (e²·a) / (2·r)
   - Итерация: b = c + s₁, s₂ = asin(p₁·sin(2b) / ✓(1 - e²·sin²b))
4. Вычисление высоты: h = p·cos φ + z·sin φ - N·(1 - e²·sin² φ)

### 2. Преобразования между датумами (по ГОСТ 32453-2017)

#### 2.1 ПЗ-90.11 ↔ WGS-84(G1150)

```csharp
// ПЗ-90.11 → WGS-84
public static (double X, double Y, double Z) ConvertPZ90ToWGS84(
    (double X, double Y, double Z) pz90)

// WGS-84 → ПЗ-90.11  
public static (double X, double Y, double Z) ConvertWGS84ToPZ90(
    (double X, double Y, double Z) wgs84)
```

#### 2.2 ПЗ-90.11 ↔ ГСК-2011

```csharp
// ПЗ-90.11 → ГСК-2011
public static (double X, double Y, double Z) ConvertPZ90ToGSK2011(
    (double X, double Y, double Z) pz90)

// ГСК-2011 → ПЗ-90.11
public static (double X, double Y, double Z) ConvertGSK2011ToPZ90(
    (double X, double Y, double Z) gsk2011)
```

#### 2.3 ПЗ-90.11 ↔ СК-42/95

```csharp
// ПЗ-90.11 → СК
public static (double X, double Y, double Z) ConvertPZ90ToSK(
    (double X, double Y, double Z) pz90)

// СК → ПЗ-90.11
public static (double X, double Y, double Z) ConvertSKToPZ90(
    (double X, double Y, double Z) sk)
```

### 3. Проекция Гаусса-Крюгера (эллипсоид Красовского)

#### 3.1 Из геодезических WGS-84 в плоские координаты Гаусса-Крюгера

```csharp
public static (double x, double y) ConvertToGaussKrueger(
    double latitudeWgs84,   // широта в радианах (WGS-84)
    double longitudeWgs84)  // долгота в радианах (WGS-84)
```

**Цепочка преобразований:**
1. WGS-84 (геодезические) → Декартовы WGS-84
2. Декартовы WGS-84 → Декартовы ПЗ-90.11
3. Декартовы ПЗ-90.11 → Декартовы СК
4. Декартовы СК → Геодезические СК (эллипсоид Красовского)
5. Проекция Гаусса-Крюгера (формулы 24-26 ГОСТ)

**Результат:**
- `x` — северное смещение (ордината) в метрах
- `y` — восточное смещение с номером зоны (первые цифры — номер зоны)

#### 3.2 Из плоских координат Гаусса-Крюгера в геодезические WGS-84

```csharp
public static (double B, double L, double H) ConvertFromGaussKrueger(
    double x,  // северное смещение в метрах
    double y)  // восточное смещение с номером зоны
```

**Обратная цепочка преобразований:**
1. Обратная проекция Гаусса-Крюгера (формулы 29-36 ГОСТ) → Геодезические СК
2. Геодезические СК → Декартовы СК
3. Декартовы СК → Декартовы ПЗ-90.11
4. Декартовы ПЗ-90.11 → Декартовы WGS-84
5. Декартовы WGS-84 → Геодезические WGS-84

### 4. Координаты локальной касательной плоскости (ENU)

#### 4.1 Из геодезических WGS-84 в локальные ENU

```csharp
// Преобразование по отдельным параметрам
public static (double east, double north, double up) ConvertGeodeticToLocalTangent(
    double targetLatitude, double targetLongitude, double targetAltitude,  // целевая точка
    double originLatitude, double originLongitude, double originAltitude)  // опорная точка

// Преобразование с использованием кортежей
public static (double east, double north, double up) ConvertGeodeticToLocalTangent(
    (double B, double L, double H) target,
    (double B, double L, double H) origin)
```

**Матрица преобразования (ECEF → ENU):**
```
E = -sin λ·dx + cos λ·dy
N = -sin φ·cos λ·dx - sin φ·sin λ·dy + cos φ·dz
U =  cos φ·cos λ·dx + cos φ·sin λ·dy + sin φ·dz
```
где (dx, dy, dz) — разность координат в ECEF, (φ, λ) — координаты опорной точки.

#### 4.2 Из локальных ENU в геодезические WGS-84

```csharp
// Преобразование по отдельным параметрам
public static (double B, double L, double H) ConvertLocalTangentToGeodetic(
    double east, double north, double up,  // локальные координаты
    double originLatitude, double originLongitude, double originAltitude)  // опорная точка

// Преобразование с использованием кортежей
public static (double B, double L, double H) ConvertLocalTangentToGeodetic(
    (double east, double north, double up) local,
    (double B, double L, double H) origin)
```

**Обратное преобразование (ENU → ECEF):**
```
dx = -sin λ·E - sin φ·cos λ·N + cos φ·cos λ·U
dy =  cos λ·E - sin φ·sin λ·N + cos φ·sin λ·U
dz =  cos φ·N + sin φ·U
```

### Особенности реализации

1. **Все угловые величины работают в радианах** — требуется предварительное преобразование, если исходные данные в градусах
2. **Итеративный алгоритм** для обратного преобразования обеспечивает высокую точность (ε = 1e-10)
3. **Цепочки преобразований** используют промежуточную систему ПЗ-90.11 как базовую для всех трансформаций
4. **Проекция Гаусса-Крюгера** реализована строго по формулам ГОСТ с эмпирическими коэффициентами для эллипсоида Красовского

### Пример использования

```csharp
// Преобразование геодезических координат WGS-84 в проекцию Гаусса-Крюгера
double lat = 55.7558 * Math.PI / 180.0;  // широта Москвы в радианах
double lon = 37.6176 * Math.PI / 180.0;  // долгота Москвы в радианах

var (x, y) = TranslationWGS84.ConvertToGaussKrueger(lat, lon);
Console.WriteLine($"X: {x:F2} м, Y: {y:F2} м");

// Вычисление расстояния между двумя точками
var point1 = (55.7558 * Math.PI / 180.0, 37.6176 * Math.PI / 180.0, 150.0);
var point2 = (55.7614 * Math.PI / 180.0, 37.6289 * Math.PI / 180.0, 160.0);

double distance = TranslationWGS84.CalculateDistanceInLocalTangent(point1, point2);
Console.WriteLine($"Расстояние: {distance:F2} м");

double azimuth = TranslationWGS84.CalculateAzimuthInLocalTangent(point1, point2);
Console.WriteLine($"Азимут: {azimuth * 180.0 / Math.PI:F2}°");
```
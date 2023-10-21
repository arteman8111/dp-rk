// ИД для всех
const m0 = 3015 // кг ~ начальная масса СВ
const m_pg = 1550 // кг ~ Предельная масса рабочего запаса топлива
const m_const = 615 // кг ~ Масса конструкции
const tv = 14 // c ~ продолжительность вертикального участка выведения
const step = 0.1 // с ~ шаг инт
const uM = 4903 * math.pow(10,9) // м^3/с^2 ~ гравитационная постоянная Луны
const rM = 1738 * math.pow(10,3) // м ~ радиус Луны

// ИД мои
const W = 3.360 * math.pow(10,3) // м/с ~ Эффективная скорость истечения газов
const P1 = 8.30 * math.pow(10,3) // Н ~ Тяга СВ 1 АУТ
const h_isl_1_1 = 140 * math.pow(10,3)// м ~ высота орбиты для P1
const h_isl_1_2 = 180 * math.pow(10,3)// м ~ высота орбиты для P1
const P2 = 9.90 * math.pow(10,3) // Н ~ Тяга СВ 2 АУТ
const h_isl_2_1 = 180 * math.pow(10,3)// м ~ высота орбиты для P2
const h_isl_2_2 = 310 * math.pow(10,3) // м ~ высота орбиты для P2

// Подгоняемые параметры
const t1 = 370.56; // с
const t2 = 450.11; // с
const thet_torch = -0.001; // рад
const thet_2 = -0.3;

// Ошибка по ТЗ
const eps_r = math.pow(10, -8)
const eps_v = math.pow(10, -8)
const eps_thet = math.pow(10, -5) * math.pi / 180;
const eps_extr = math.pow(10, -2)

// НУ
const t0 = 0;
const vx0 = 0;
const vy0 = 0;
const x0 = 0;
const y0 = 0;
const v0 = 0;
const r0 = rM;
const thet = math.pi / 2;
const THET = 0;
const THETc = math.pi / 2;
const alfa = 0;
const fi = 0;

export {
    m0,
    m_pg,
    m_const,
    tv,
    step,
    uM,
    rM,
    W,
    P1,
    P2,
    h_isl_1_1,
    h_isl_1_2,
    h_isl_2_1,
    h_isl_2_2,
    t1,
    t2,
    thet_torch,
    eps_r,
    eps_v,
    eps_thet,
    eps_extr,
    thet_2,
    t0,
    vx0,
    vy0,
    x0,
    y0,
    v0,
    r0,
    THET,
    THETc,
    alfa,
    fi,
    thet
}
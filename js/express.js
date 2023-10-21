import * as param from "./const.js";

const thet = (t, thet_torch, thet_2, t1, t2) => {
    if (t <= param.tv) {
        return math.pi / 2
    } else if (t <= t1) {
        return math.pi / 2 + thet_torch * (t - param.tv)
    } else if (t < t2) {
        return math.pi / 2 + thet_torch * (t1 - param.tv)
    } else if (t >= t2) {
        return thet_2
    }
}
const P = (t, t1, t2, P0) => {
    if (t <= t1) {
        return P0;
    } else if (t <= t2 && t > t1) {
        return 0;
    } else if (t > t2) {
        return P0;
    }
}
const TETA = (vx, vy, x, y, v, r) => math.asin((x * vx + (y + param.rM) * vy) / (v * r));
const TETAc = (vx, vy) => math.atan(vy / vx);
const alpha = (thet, THET) => thet - THET;
const fi = (x, y) => math.atan(x / (y + param.rM));

const gx = (x, y) => param.uM * x / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)
const gy = (x, y) => param.uM * (param.rM + y) / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)

const vx_toch = (t, m, x, y, P, thet_torch, thet_2, t1, t2) => P * math.cos(thet(t, thet_torch, thet_2, t1, t2)) / m - gx(x, y);
const vy_toch = (t, m, x, y, P, thet_torch, thet_2, t1, t2) => P * math.sin(thet(t, thet_torch, thet_2, t1, t2)) / m - gy(x, y);
const x_toch = (v) => v;
const y_toch = (v) => v;
const m_toch = (P) => -P / param.W;

// Конечные параметры
const vk = (h) => math.sqrt(param.uM / (param.rM + h));
const rk = (h) => h + param.rM;

export {
    vx_toch,
    vy_toch,
    x_toch,
    y_toch,
    m_toch,
    P,
    thet,
    TETA,
    TETAc,
    alpha,
    fi,
    vk,
    rk
}
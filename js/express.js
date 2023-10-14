import * as param from "./const.js";

const thet = (t) => {
    if (t <= param.tv) {
        return math.pi / 2
    } else if (t <= param.t1) {
        return math.pi / 2 + param.thet_torch * (t - param.tv)
    } else if (t < param.t2) {
        return math.pi / 2 + param.thet_torch * (param.t1 - param.tv)
    } else if (t >= param.t2) {
        return param.thet_2
    }
}
const P = (t) => {
    if (t <= param.t1){
        return param.P2;
    } else if (t <= param.t2 && t > param.t1) {
        return 0;
    } else if (t > param.t2){
        return param.P2;
    }
}
const TETA = (vx, vy, x, y, v, r) => math.asin((x * vx + y * vy) / (v * r));
const TETAc = (vx, vy) => math.atan(vy / vx);
const alpha = (vy,y) => math.acos(y/vy);
const fi = (x, y) => math.atan(x / (y + param.rM));

const gx = (x, y) => param.uM * x / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)
const gy = (x, y) => param.uM * (param.rM + y) / math.pow((math.pow(x, 2) + math.pow((param.rM + y), 2)), 1.5)

const vx_toch = (t, m, x, y, P) => P * math.cos(thet(t)) / m - gx(x, y);
const vy_toch = (t, m, x, y, P) => P * math.sin(thet(t)) / m - gy(x, y);
const x_toch = (v) => v;
const y_toch = (v) => v;
const m_toch = (P) => -P / param.W

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
    fi
}
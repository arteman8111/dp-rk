import { rM } from "./const.js";
import { vx_toch, vy_toch, x_toch, y_toch, m_toch} from "./express.js";
import { P, thet, TETA, TETAc, alpha, fi } from "./express.js";
const grad = value => value * 180 / math.pi;
function rungekutta(el, dt, t) {
    // debugThis(t, param.t1)   
    let k1_vx = dt * vx_toch(t, el[0], el[3], el[4], P(t));
    let k1_vy = dt * vy_toch(t, el[0], el[3], el[4], P(t));
    let k1_x = dt * x_toch(el[1]);
    let k1_y = dt * y_toch(el[2]);
    let k1_m = dt * m_toch(P(t))

    let k2_vx = dt * vx_toch(t + dt / 2, el[0] + k1_m / 2, el[3] + k1_x / 2, el[4] + k1_y / 2, P(t));
    let k2_vy = dt * vy_toch(t + dt / 2, el[0] + k1_m / 2, el[3] + k1_x / 2, el[4] + k1_y / 2, P(t));
    let k2_x = dt * x_toch(el[1] + k1_vx / 2);
    let k2_y = dt * y_toch(el[2] + k1_vy / 2);
    let k2_m = dt * m_toch(P(t));

    let k3_vx = dt * vx_toch(t + dt / 2, el[0] + k2_m / 2, el[3] + k2_x / 2, el[4] + k2_y / 2, P(t));
    let k3_vy = dt * vy_toch(t + dt / 2, el[0] + k2_m / 2, el[3] + k2_x / 2, el[4] + k2_y / 2, P(t));
    let k3_x = dt * x_toch(el[1] + k2_vx / 2);
    let k3_y = dt * y_toch(el[2] + k2_vy / 2);
    let k3_m = dt * m_toch(P(t));

    let k4_vx = dt * vx_toch(t + dt, el[0] + k3_m, el[3] + k3_x, el[4] + k3_y, P(t));
    let k4_vy = dt * vy_toch(t + dt, el[0] + k3_m, el[3] + k3_x, el[4] + k3_y, P(t));
    let k4_x = dt * x_toch(el[1] + k3_vx);
    let k4_y = dt * y_toch(el[2] + k3_vy);
    let k4_m = dt * m_toch(P(t))


    el[0] += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6;
    el[1] += (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx) / 6;
    el[2] += (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy) / 6;
    el[3] += (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    el[4] += (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    el[5] = math.sqrt(math.pow(el[1], 2) + math.pow(el[2], 2))
    el[6] = math.sqrt(math.pow(el[3], 2) + math.pow(el[4] + rM, 2));
    const TETAk = TETA(el[1], el[2], el[3], el[4], el[5], el[6]);
    const TETAkc = TETAc(el[1], el[2])
    el[7] = grad(thet(t));  
    el[8] = grad(TETAkc);
    el[9] = grad(alpha(TETAkc, TETAk));
    el[10] = grad(fi(el[3], el[4]));
}
export default rungekutta
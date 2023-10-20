import * as param from "./const.js"
import rungekutta from "./rungekutta.js";
import { html } from "./render.js";

const integr = (thet) => {
    // html(t, el)
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0, param.thet, param.THET, param.alfa, param.fi, param.THETc];
    
    let t_prev = t;
    let el_prev = el.slice();
    while (math.abs(el[5] - param.vk) > param.eps_v) {
        if (t + dt > param.t1 && t < param.t1) {
            dt = param.t1 - t;
            t += dt;
            rungekutta(el, dt, t, thet[0], thet[1]);
            // html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t, thet[0], thet[1]);
            // html(t, el);
            dt = param.step
        }
        if (t + dt > param.t2 && t < param.t2) {
            dt = param.t2 - t;
            t += dt;
            rungekutta(el, dt, t, thet[0], thet[1]);
            // html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t, thet[0], thet[1]);
            // html(t, el)
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        rungekutta(el, dt, t, thet[0], thet[1]);
        t += dt;
        // html(t, el);
        if (el[5] > param.vk) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    return el
}

function init() {
    let thet = [param.thet_torch, param.thet_2];
    const thet_step = [math.pow(10, -5), math.pow(10, -8)];
    
    let el1, delta_F;

    function thet_torch_rad(thet) {
        let el, thetk, delta_next, delta_prev;
        thetk = thet.slice();
        thetk[0] += thet_step[0];
        el = integr(thetk);
        delta_next = el[6] - param.rk;
        thetk[0] -= 2 * thet_step[0];
        el = integr(thetk);
        delta_prev = el[6] - param.rk
        return (delta_next - delta_prev) / (2 * thet_step[0])
    }
    function thet_torch_th(thet) {
        let el, thetk, delta_next, delta_prev;
        thetk = thet.slice();
        thetk[0] += thet_step[0];
        el = integr(thetk);
        delta_next = el[8];
        thetk[0] -= 2 * thet_step[0];
        el = integr(thetk);
        delta_prev = el[8];
        return (delta_next - delta_prev) / (2 * thet_step[0])
    }
    function thet_2_rad(thet) {
        let el, thetk, delta_next, delta_prev;
        thetk = thet.slice();
        thetk[1] += thet_step[1];
        el = integr(thetk);
        delta_next = el[6] - param.rk;
        thetk[1] -= 2 * thet_step[1];
        el = integr(thetk);
        delta_prev = el[6] - param.rk;
        return (delta_next - delta_prev) / (2 * thet_step[1])
    }
    function thet_2_th(thet) {
        let el, thetk, delta_next, delta_prev;
        thetk = thet.slice();
        thetk[1] += thet_step[1];
        el = integr(thetk);
        delta_next = el[8];
        thetk[1] -= 2 * thet_step[1];
        el = integr(thetk);
        delta_prev = el[8];
        return (delta_next - delta_prev) / (2 * thet_step[1])
    }
     do {
         el1 = integr(thet);
         delta_F = [-el1[6] + param.rk, -el1[8]];
         if(math.sqrt(math.pow((delta_F[0])/param.eps_r,2) + math.pow((delta_F[1])/param.eps_thet,2)) < 1){
            break
         }
         const J = math.inv(math.matrix([[thet_torch_rad(thet), thet_2_rad(thet)],[thet_torch_th(thet), thet_2_th(thet)]]));
         const U = math.multiply(J, delta_F);
         thet = math.add(thet, U)._data;
     } while (true)
     const el2 = integr(thet);
     console.log(el2[5] - param.vk, el2[6] - param.rk, el2[8]);
}
init()
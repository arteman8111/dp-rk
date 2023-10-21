import * as param from "./const.js"
import paramIter from "./paramIter.js";
import { html, render } from "./render.js";
import { vk, rk } from "./express.js";
let i = 1;
const integr = (thet, h, t1, t2, P0) => {
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0, param.thet, param.THET, param.alfa, param.fi, param.THETc];

    let t_prev = t;
    let el_prev = el.slice();
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t2) {
            dt = t1 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            dt = param.step;
            t += dt;
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            dt = param.step;
            t += dt;
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
        t += dt;
        if (el[5] > vk(h)) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    return el
}
const printTable = (thet, h, t1, t2, P0) => {
    render(thet, P0, h, t1, t2, i);
    let dt = param.step;
    let t = param.t0;
    let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0, param.thet, param.THET, param.alfa, param.fi, param.THETc];
    html(t, el, i);

    let t_prev = t;
    let el_prev = el.slice();
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t1) {
            dt = t1 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i);
            dt = param.step
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step - dt;
            t += dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
            html(t, el, i)
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0);
        t += dt;
        if (el[5] < vk(h) && dt === param.step) {
            html(t, el, i);
        }
        if (math.abs(el[5] - vk(h)) < param.eps_v) {
            html(t, el, i);
        }
        if (el[5] > vk(h)) {
            el = el_prev.slice();
            t = t_prev;
            el_prev = el.slice();
            t_prev = t;
            dt = dt / 10;
        }
    }
    i++
}


function init() {
    function optimus() {
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
            if (math.sqrt(math.pow((delta_F[0]) / param.eps_r, 2) + math.pow((delta_F[1]) / param.eps_thet, 2)) < 1) {
                break
            }
            const J = math.inv(math.matrix([[thet_torch_rad(thet), thet_2_rad(thet)], [thet_torch_th(thet), thet_2_th(thet)]]));
            const U = math.multiply(J, delta_F);
            thet = math.add(thet, U)._data;
        } while (true)
        return thet
    }
    const thet_id = [param.thet_torch, param.thet_2];
    const el3 = printTable(thet_id, param.h_isl_2_2, param.t1, param.t2, param.P2);
}
init()
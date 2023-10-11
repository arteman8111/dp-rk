import * as param from "./const.js"
import rungekutta from "./rungekutta.js";
import { html } from "./render.js";

// НУ
let dt = param.step;
let t = param.t0;
let t_prev = t;
let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0,0,0,0,0];
let el_prev = el.slice();

function init() {
    html(t, el)
    while (
        el[5] < param.vk
        // t < 20
        ) {
        if (t + dt > param.t1 && t < param.t1) {
            dt = param.t1 - t;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el);
            dt = param.step
        }
        if (t + dt > param.t2 && t < param.t2) {
            dt = param.t2 - t;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        rungekutta(el, dt, t);
        t += dt;
        if (el[5] < param.vk) {
            html(t, el);
        }
        if (el[5] - param.eps_v > param.vk) {
            el = el_prev.slice();
            t = t_prev;
            dt = dt / 10;
            t += dt;
            rungekutta(el, dt, t);
        }
    }
}
init()
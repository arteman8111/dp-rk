import * as param from "./const.js"
import rungekutta from "./rungekutta.js";
import { html } from "./render.js";
import { secantMethod } from "./secantMethod.js";

// НУ
let dt = param.step;
let t = param.t0;
let t_prev = t;
let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0, 90, 0, 0, 0, 90];
let el_prev = el.slice();

function init() {
    html(t, el)
    while (el[5] + param.eps_v < param.vk) {
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
        if (el[5] + param.eps_v > param.vk) {
            el = el_prev.slice();
            el_prev = el.slice();
            t = t_prev;
            t_prev = t;
            dt = dt / 10;
            rungekutta(el, dt, t);
            t += dt;
            if (el[5] + param.eps_v > param.vk) {
                el = el_prev.slice();
                t = t_prev;
            }
        }
        // if(dt === param.step){
        //     html(t, el);
        // } else if (math.abs(el[5] - param.vk) <= param.eps_v){
        //     html(t,el)
        // }
    }
}
init()
console.log(el);
console.log(el[5] - param.vk);
const delta_r = el[6] - param.rk;
const delta_thet = el[8]
let err_arr = [delta_r,delta_thet]
console.log(err_arr);

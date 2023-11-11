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
    let boolk = false;
    while (math.abs(el[5] - vk(h)) > param.eps_v) {
        if (t + dt > t1 && t < t1) {
            dt = t1 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            boolk = true;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            dt = param.step;
        }
        if (t + dt > t2 && t < t2) {
            dt = t2 - t;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            boolk = false;
            dt = param.step - dt;
            paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
            t += dt;
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        paramIter(el, dt, t, thet[0], thet[1], t1, t2, P0, boolk);
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
const thet_optimus = (thet_arr, t, P0, h) => {
    let thet = thet_arr.slice();
    const thet_step = [math.pow(10, -5), math.pow(10, -8)];
    let el1, dF;
    function thet_1_2_torch(thet) {
        let el, thetk, delta_next_r_1, delta_next_r_2, delta_next_thet_1, delta_next_thet_2;
        let delta_prev_r_1, delta_prev_r_2, delta_prev_thet_1, delta_prev_thet_2;
        thetk = thet.slice();
        thetk[0] += thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_1 = el[6] - rk(h);
        delta_next_thet_1 = el[8];

        thetk[0] -= 2 * thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_1 = el[6] - rk(h);
        delta_prev_thet_1 = el[8];

        thetk = thet.slice();
        thetk[1] += thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_2 = el[6] - rk(h);
        delta_next_thet_2 = el[8];

        thetk[1] -= 2 * thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_2 = el[6] - rk(h);
        delta_prev_thet_2 = el[8];

        const dr_1 = (delta_next_r_1 - delta_prev_r_1) / (2 * thet_step[0]);
        const dthet_1 = (delta_next_thet_1 - delta_prev_thet_1) / (2 * thet_step[0]);
        const dr_2 = (delta_next_r_2 - delta_prev_r_2) / (2 * thet_step[1]);
        const dthet_2 = (delta_next_thet_2 - delta_prev_thet_2) / (2 * thet_step[1]);
        return [[dr_1, dr_2], [dthet_1, dthet_2]]
    }
    do {
        el1 = integr(thet, h, t[0], t[1], P0);
        dF = [-el1[6] + rk(h), -el1[8]];
        if (math.sqrt(math.pow((dF[0]) / param.eps_r, 2) + math.pow((dF[1]) / param.eps_thet, 2)) < 1) {
            break
        }
        const J = math.inv(math.matrix(thet_1_2_torch(thet)));
        const U = math.multiply(J, dF);
        thet = math.add(thet, U)._data;
    } while (true)
    return thet
}
const optimus = (thet_arr, t, P0, h) => {
    let thet = thet_arr.slice();
    const thet_step = [math.pow(10, -5), math.pow(10, -8)];
    let el1, dF;
    function thet_1_2_torch(thet) {
        let el, thetk, delta_next_r_1, delta_next_r_2, delta_next_thet_1, delta_next_thet_2;
        let delta_prev_r_1, delta_prev_r_2, delta_prev_thet_1, delta_prev_thet_2;
        thetk = thet.slice();
        thetk[0] += thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_1 = el[6] - rk(h);
        delta_next_thet_1 = el[8];

        thetk[0] -= 2 * thet_step[0];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_1 = el[6] - rk(h);
        delta_prev_thet_1 = el[8];

        thetk = thet.slice();
        thetk[1] += thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_next_r_2 = el[6] - rk(h);
        delta_next_thet_2 = el[8];

        thetk[1] -= 2 * thet_step[1];
        el = integr(thetk, h, t[0], t[1], P0);
        delta_prev_r_2 = el[6] - rk(h);
        delta_prev_thet_2 = el[8];

        const dr_1 = (delta_next_r_1 - delta_prev_r_1) / (2 * thet_step[0]);
        const dthet_1 = (delta_next_thet_1 - delta_prev_thet_1) / (2 * thet_step[0]);
        const dr_2 = (delta_next_r_2 - delta_prev_r_2) / (2 * thet_step[1]);
        const dthet_2 = (delta_next_thet_2 - delta_prev_thet_2) / (2 * thet_step[1]);
        return [[dr_1, dr_2], [dthet_1, dthet_2]]
    }
    do {
        el1 = integr(thet, h, t[0], t[1], P0);
        dF = [-el1[6] + rk(h), -el1[8]];
        if (math.sqrt(math.pow((dF[0]) / param.eps_r, 2) + math.pow((dF[1]) / param.eps_thet, 2)) < 1) {
            break
        }
        const J = math.inv(math.matrix(thet_1_2_torch(thet)));
        const U = math.multiply(J, dF);
        thet = math.add(thet, U)._data;
    } while (true)
    return el1[0]
}
function init() {
    function levenberg() {
        function dm_t1(t0_arr) {
            let delta_next, delta_prev, tk_arr;
            tk_arr = t0_arr.slice();
            tk_arr[0] += dt_step;
            delta_next = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            tk_arr[0] -= 2 * dt_step;
            delta_prev = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            return (delta_next - delta_prev) / (2 * dt_step)
        }

        function dm_t2(t0_arr) {
            let delta_next, delta_prev, tk_arr;
            tk_arr = t0_arr.slice();
            tk_arr[1] += dt_step;
            delta_next = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            tk_arr[1] -= 2 * dt_step;
            delta_prev = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            return (delta_next - delta_prev) / (2 * dt_step)
        }

        function ddm_t1(t0_arr) {
            let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
            tk_arr = t0_arr.slice();
            dm1_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            dm2_1 = dm1_2;
            tk_arr[0] += 2 * dt_step_double;
            dm2_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            tk_arr[0] -= 4 * dt_step_double;
            dm1_1 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            z1 = (dm1_1 - dm1_2) / (2 * dt_step_double);
            z2 = (dm2_1 - dm2_2) / (2 * dt_step_double);
            f_res = (z1 - z2) / (2 * dt_step_double);
            return f_res
        }

        function ddm_t2(t0_arr) {
            let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
            tk_arr = t0_arr.slice();
            dm1_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            dm2_1 = dm1_2;
            tk_arr[1] += 2 * dt_step_double;
            dm2_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            tk_arr[1] -= 4 * dt_step_double;
            dm1_1 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            z1 = (dm1_1 - dm1_2) / (2 * dt_step_double);
            z2 = (dm2_1 - dm2_2) / (2 * dt_step_double);
            f_res = (z1 - z2) / (2 * dt_step_double);
            return f_res
        }
        function ddm_t1_t2(t0_arr) {
            let z1, z2, f_res, dm1_1, dm1_2, dm2_1, dm2_2, tk_arr;
            tk_arr = t0_arr.slice();
            // tk_arr.map(el => el - dt_step_double);
            tk_arr[0] += dt_step;
            tk_arr[1] += dt_step_double;
            dm1_1 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            // tk_arr.map(el => el + 2 * dt_step_double);
            tk_arr[0] -= 2 * dt_step;
            dm2_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            // tk_arr[1] -= 2 * dt_step_double;
            tk_arr[0] += 2 * dt_step;
            tk_arr[1] -= 2 * dt_step_double;
            dm2_1 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);
            tk_arr[0] -= 2 * dt_step;
            dm1_2 = optimus(thet_id, tk_arr, param.P2, param.h_isl_2_2);

            z1 = (dm1_1 - dm2_2) / (2 * dt_step);
            z2 = (dm2_1 - dm1_2) / (2 * dt_step);
            f_res = (z1 - z2) / (2 * dt_step_double);
            return f_res
        }

        let thet_id;
        let t_id = [param.t1, param.t2];
        let dt_step = math.pow(10, -2);
        let dt_step_double = 2 * math.pow(10, -2);
        
        const I = math.matrix(
            [
                [1, 0],
                [0, 1]
            ]
        )
        let C1 = math.pow(2, -1);
        let C2 = math.pow(2, 1);
        let alfa = math.pow(2, 1);
        let alfa_iter = alfa;
        let t_iter = t_id.slice();
        let t_prev = t_iter.slice();

        let gradient, gradient_module
        let m_iter, m_iter_prev
        thet_id = thet_optimus([param.thet_torch, param.thet_2], t_iter, param.P2, param.h_isl_2_2);
        do {
            gradient = math.matrix([dm_t1(t_iter), dm_t2(t_iter)]);
            gradient_module = math.sqrt(math.pow(gradient._data[0], 2) + math.pow(gradient._data[1], 2));
            m_iter = optimus(thet_id, t_iter, param.P2, param.h_isl_2_2);
            console.log("Время: ", t_iter);
            console.log("Градиент: ", gradient_module);
            console.log("Mass", m_iter);
            console.log("alfa", alfa_iter);
            if (gradient_module < param.eps_extr) {
                break
            }
            // Гессиан в цикле пока модуль больше предыдушего значения модуля
            // gradient_prev = math.clone(gradient);
            m_iter_prev = m_iter;
            t_prev = t_iter.slice();
            let gessian = math.matrix(
                [
                    [ddm_t1(t_iter), ddm_t1_t2(t_iter)],
                    [ddm_t1_t2(t_iter), ddm_t2(t_iter)]
                ]
            )
            console.log('gessian',gessian);
            do {
                let first = math.add(gessian, math.multiply(alfa_iter, I)); // (Hi + alfai * I)
                let second = math.inv(first); // // (Hi + alfai * I)^-1
                console.log('second',second);
                let third = math.multiply(second, gradient); // (Hi + alfai * I)^-1 x gradient
                console.log('third', third);
                // third = math.multiply(third, -1);
                t_iter = math.add(t_iter, third)._data;
                // gradient = math.matrix([dm_t1(t_iter), dm_t2(t_iter)]);
                // gradient_module = math.sqrt(math.pow(gradient._data[0], 2) + math.pow(gradient._data[1], 2));
                m_iter = optimus(thet_id, t_iter, param.P2, param.h_isl_2_2);

                if (m_iter > m_iter_prev) {
                    alfa_iter *= C1;
                    break;
                }
                alfa_iter *= C2;
                t_iter = t_prev.slice();
                // gradient = math.clone(gradient_prev);
            } while (true);
        } while (true);
    }
    levenberg()
}
init()
